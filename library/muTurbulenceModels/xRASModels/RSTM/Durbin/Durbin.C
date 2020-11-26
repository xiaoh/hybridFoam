/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "Durbin.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"
#include "fixedInternalValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Durbin, 0);
addToRunTimeSelectionTable(RASModel, Durbin, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Durbin::Durbin
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),
    GenElliptic(U, phi, lamTransportModel),

    solveK_(coeffDict_.lookupOrAddDefault<Switch>("solveK", true)),
    fBC_(coeffDict_.lookupOrAddDefault<word>("fBC", "automatic")),
    durFlag_(1.0),
    crossTurbDiffusion_(coeffDict_.lookupOrAddDefault<Switch>("crossTurbDiffusion", false)),
    zeroTraceF_(coeffDict_.lookupOrAddDefault<Switch>("zeroTraceF", false)),
    wallsAlignedWithZ_(coeffDict_.lookupOrAddDefault<Switch>("wallsAlignedWithZ", true)),

    ellipticOperatorCorrection_(coeffDict_.lookupOrAddDefault<Switch>("ellipticOperatorCorrection", false)),
    beta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta",
            coeffDict_,
            0.083
        )
    ),

    f_
    (
        IOobject
        (
            "f",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{
    if(fBC_ == "Hanjalic")
    {
        durFlag_ = 0.0;
    }
    else if(fBC_ == "Durbin")
    {
        durFlag_ = 1.0;
    }
    else // "automatic"
    {
        durFlag_ = scalar(int(bool(solveK_)));
    }
    coeffDict_.set<scalar>("durFlag", durFlag_);
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Durbin::read()
{
    if (GenElliptic::read())
    {
        solveK_.readIfPresent(word("solveK"), coeffDict());
        ellipticOperatorCorrection_.readIfPresent(word("ellipticOperatorCorrection"), coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


void Durbin::correct()
{
    GenElliptic::correct();

    if (!turbulence_)
    {
        return;
    }

    volSymmTensorField P = -twoSymm(R_ & fvc::grad(U_));
    volScalarField G("RASModel::G", 0.5*mag(tr(P)));

    volScalarField Ts("T", T());

    #include "../include/epsilonWallI2.H" // set patch internal eps values

    // split R_ into normal diffusion and cross diffusion terms
    volSymmTensorField Rdiag = R_;
    dimensionedScalar kzero = k0_ * 0.0;
    Rdiag.replace(symmTensor::XY, kzero);
    Rdiag.replace(symmTensor::YZ, kzero);
    Rdiag.replace(symmTensor::XZ, kzero);
    volSymmTensorField Rupper = R_ - Rdiag;

    symmTensor  minDiagR = gMin(Rdiag);
    if(
        minDiagR.xx() < 0.0 ||
        minDiagR.yy() < 0.0 ||
        minDiagR.zz() < 0.0 )
    {
        Info << "muDurbin::correct():: Warning! " << nl
             << "negative diagonal for R. I will probably fail soon! Rdiag.min = "
             << minDiagR << endl;
    }

    if(debug)
    {
        Info << "  max(C2/T): " << gMax((C2_/Ts)()) << endl;
    }

    surfaceScalarField Tsf = fvc::interpolate(Ts, "interpolate(T)");
    surfaceSymmTensorField Rdiagf  = fvc::interpolate(Rdiag, "interpolate(R)");
    surfaceSymmTensorField Rupperf = fvc::interpolate(Rupper, "interpolate(R)");

    // Dissipation equation 
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
        + fvm::div(phi_, epsilon_)
        - fvm::Sp(fvc::div(phi_), epsilon_)
        - fvm::laplacian(Cmu_/sigmaEps_ * Tsf * Rdiagf, epsilon_, "laplacian(epsilon)")
        - fvm::laplacian(nu(), epsilon_, "laplacian(epsilon)")
      ==
        C1_ * G/Ts * ( 1.0 + 0.1*G/epsilon_)
        - fvm::Sp(C2_/Ts, epsilon_)
    );

    if(crossTurbDiffusion_)
    {
        epsEqn() -= fvc::laplacian(Cmu_/sigmaEps_ * Tsf * Rupperf, epsilon_, "laplacian(epsilon)");
    }

    epsEqn().relax();
    epsEqn().boundaryManipulate(epsilon_.boundaryField());
    solve(epsEqn);
    bound(epsilon_, epsilon0_);

    // TKE equation
    if(solveK_)
    {
        tmp<fvScalarMatrix> kEqn
            (
                fvm::ddt(k_)
                + fvm::div(phi_, k_)
                - fvm::Sp(fvc::div(phi_), k_)
                - fvm::laplacian(Cmu_/sigmaK_ * Tsf * Rdiagf, k_, "laplacian(k)")
                - fvm::laplacian(nu(), k_, "laplacian(k)")
                ==
                G
                - fvm::Sp(epsilon_/k_, k_)
            );

        if(crossTurbDiffusion_)
        {
            kEqn() -= fvc::laplacian(Cmu_/sigmaK_ * Tsf * Rupperf, k_, "laplacian(k)");
        }
        
        kEqn().relax();
        solve(kEqn);
    }
    else
    {
        k_ = 0.5 * tr(R_);
    }
    bound(k_, k0_);

    // Reynolds stress equation
    #include "fWallI.H" // set patch internal f values

    if(debug)
    {
        Info << "  max(divPhi): " << gMax((fvc::div(phi_))()) << endl;
        Info << "  max(epsilon/k) : " << gMin((epsilon_/k_)()) << endl;
        Info << "  min(k*f): " << gMin((k_ * f_)()) << endl;
        Info << "  min tr(k*f): " << gMin(tr(k_ * f_)()) << endl;
        Info << "  min tr(f): " << gMin(tr(f_)()) << endl;
        Info << "  min(P + k*f): " << gMin((k_ * f_ + P)()) << endl;
        Info << "  min tr(P + k*f): " << gMin(tr(k_ * f_ + P)()) << endl;
        Info << "  min(diagR): (" << minDiagR.xx() << ", " << minDiagR.yy() 
             << ", " << minDiagR.zz() << ")" << endl;
    }
    
    tmp<fvSymmTensorMatrix> REqn
        (
            fvm::ddt(R_)
            + fvm::div(phi_, R_)
            - fvm::Sp(fvc::div(phi_), R_)
            - fvm::laplacian(Cmu_/sigmaK_ * Tsf * Rdiagf, R_, "laplacian(R)")
            - fvm::laplacian(nu(), R_, "laplacian(R)")
            + fvm::Sp(epsilon_/k_, R_)
            ==                                        
            P                                         // production tensor
            + k_ * f_
        );

    if(crossTurbDiffusion_)
    {
        REqn() -= fvc::laplacian(Cmu_/sigmaK_*Ts*Rupper, R_, "laplacian(R)");
    }

    REqn().relax();
    solve(REqn);

    if(solveK_)
    {
        forAll(R_, celli)
        {
            symmTensor& rij = R_.internalField()[celli];
            rij.zz() = 2.0*k_.internalField()[celli] - rij.xx() - rij.yy();
        }
    }
    
    R_.max // bound diagonal components of R
    (
        dimensionedSymmTensor
        (
            "zero",
            R_.dimensions(),
            symmTensor
            (
                k0_.value(), -GREAT,      -GREAT,
                             k0_.value(), -GREAT,
                                          k0_.value()
            )
        )
    );

    volScalarField Ls = L();
    Ts = T(); // re-compute time scale

    volSymmTensorField exSrc = -Clrr1_*dev(R_)/Ts - Clrr2_*dev(P);

    if(SSG_)
    {
        volSymmTensorField bij("bij", 0.5 * dev(R_/k_));
        volTensorField fbij(symm2full(bij));
        volTensorField     gradU = fvc::grad(U_);
        volSymmTensorField Sij("Sij", dev(symm(gradU)));
        volTensorField fSij(symm2full(Sij));
        volTensorField     Wij =  skew(gradU);
        
        exSrc = -(Cg1_ * k_/Ts + Cg1s_ * G) * bij
            +    Cg2_ * k_/Ts * dev(symm(bij & bij))
            +   (Cg3_ - Cg3s_ * sqrt(bij && bij)) * k_ * Sij
            +    Cg4_ * k_ * dev( twoSymm(fbij & fSij) )
            +    Cg5_ * k_ * (twoSymm(bij & Wij )) ;
    }
    
    
    tmp<fvSymmTensorMatrix> fEqn
        (
            fvm::laplacian(f_)
            ==
            fvm::Sp(1.0/sqr(Ls), f_)
            -
            (
                exSrc/k_ + dev(R_)/(k_*Ts)
            ) / sqr(Ls)
            );

    if(ellipticOperatorCorrection_)
    {
        volVectorField gradL = fvc::grad(Ls);

        volSymmTensorField LdF(        
            IOobject
            (
                "LnF",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedSymmTensor("Ldf::zero", dimensionSet(0, -1, -1, 0, 0), symmTensor::zero)
        );
        
        for(label compi=0; compi<6; compi++)
        {
            volScalarField fcomp = f_.component(compi);
            LdF.replace(compi, gradL & fvc::grad(fcomp));
        }
        // Append the correction terms. (note == above is equivalent to -)
        fEqn() -= fvm::Sp((16.0*beta_*magSqr(gradL))/sqr(Ls), f_) - 8.0*beta_* LdF/Ls;
    }
    
    fEqn().relax();
    fEqn().boundaryManipulate(f_.boundaryField());
    solve(fEqn);
    
    if(zeroTraceF_)
    {
        f_ = dev(f_);
        Info <<"  After zero out, min tr(f) = " << gMin(tr(f_)()) << endl;
    }

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
