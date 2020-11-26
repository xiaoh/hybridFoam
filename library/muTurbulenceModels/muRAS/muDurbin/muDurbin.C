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

#include "muDurbin.H"
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

defineTypeNameAndDebug(muDurbin, 0);
addToRunTimeSelectionTable(muRASModel, muDurbin, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

muDurbin::muDurbin
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel,
    const volSymmTensorField & RAvg,
    const volScalarField & epsilonAvg
)
:
    muRASModel(typeName, U, phi, lamTransportModel, RAvg, epsilonAvg),
    muGenElliptic(U, phi, lamTransportModel, RAvg, epsilonAvg),

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

bool muDurbin::read()
{
    if (muGenElliptic::read())
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

//- enforce LES results on RAS fields via relax forcing
  void muDurbin::enforceFields
  (
   dimensionedScalar turbRelaxTime, scalar rampQ, word enforceMode
  )
  {
    const volScalarField& lesFlag( mesh_.lookupObject<const volScalarField>("lesFlagR") );
    volScalarField kLES = scalar(0.5) * tr(RAvg_); // Update kLES
    label timeIndex = runTime_.timeIndex();

    scalar rampQLim = Foam::min(Foam::max(rampQ, scalar(0)), scalar(1));

    volScalarField epsilonAvgBnd = max( min(epsilonAvg_, epsilon_ * boundFactor_), epsilon_ / boundFactor_);
    volSymmTensorField RAvgBnd = max( min(RAvg_, R_ * boundFactor_), R_ / boundFactor_);

    if(turbRelaxTime < runTime_.deltaT() * thresholdImpose_)
    {
        // Directly set fields
        // prepare lists for LES cell labels and values 
        // (for setting values in corresponding matrix later)
        if(solveK_)
        {
            kLES = max( min(kLES, k_ * boundFactor_), k_ / boundFactor_);
        }

        if( timeIndex % imposeTurbEvery_ == 0)
        {
            directImpose_ = true;
            Info << "Update imposed values." << endl;
            label N = sum(pos(lesFlag.internalField()-0.5));
            
            lesCellLables_.clear(); 
            lesCellValuesK_.clear();
            lesCellValuesEpsilon_.clear();
            lesCellValuesR_.clear();
            
            lesCellLables_.setSize(N); 
            lesCellValuesK_.setSize(N);
            lesCellValuesEpsilon_.setSize(N);
            lesCellValuesR_.setSize(N);
            
            label i=0;
            forAll(lesFlag, celli) // scan lesFlag Field
            {
                if(lesFlag[celli] > 0.5)
                {
                    lesCellLables_[i] = celli;
                    lesCellValuesK_[i] = kLES[celli];
                    lesCellValuesEpsilon_[i] = epsilonAvg_[celli];
                    lesCellValuesR_[i] = RAvg_[celli];
                    i++;
                }
            }
        }
        else
        {
            Info << "Skipped update imposed values: #"
                 << timeIndex % imposeTurbEvery_
                 << endl;
        }
        
        QR_ *= scalar(0.0);
        Qepsilon_ *= scalar(0.0);
      }
    else
      {
        //- Use relax forcing to indirectly influnce K/E fields
        //  through the "correction" operation of next step
        // Compute the forcing on the LES zone only
          Info << "Enforcing all turb quantities on RAS via forcing." << endl;
          QR_ = (RAvg_ - R_) / turbRelaxTime * lesFlag * rampQLim;
          Qepsilon_ = 
              (epsilonAvg_ - epsilon_) / turbRelaxTime 
              * lesFlag * rampQLim;
      }
  }


void muDurbin::correct()
{
    muGenElliptic::correct();

    if (!turbulence_)
    {
        return;
    }

    volSymmTensorField P = -twoSymm(R_ & fvc::grad(U_));
    volScalarField G("RASModel::G", 0.5*mag(tr(P)));

    volScalarField Ts("T", T());

    #include "epsilonWallI.H" // set patch internal eps values

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
    surfaceSymmTensorField Rdiagf  = fvc::interpolate(Rdiag, "interpolate(RR)");
    surfaceSymmTensorField Rupperf = fvc::interpolate(Rupper, "interpolate(RR)");

    // Dissipation equation 
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
        + fvm::div(phi_, epsilon_)
        - fvm::Sp(fvc::div(phi_), epsilon_)
        - fvm::laplacian(Cmu_/sigmaEps_ * Tsf * Rdiagf, epsilon_, "laplacian(epsilonR)")
        - fvm::laplacian(nu(), epsilon_, "laplacian(epsilonR)")
      ==
        C1_ * G/Ts * ( 1.0 + 0.1*G/epsilon_)
        - fvm::Sp(C2_/Ts, epsilon_)
        + Qepsilon_
    );

    if(crossTurbDiffusion_)
    {
        epsEqn() -= fvc::laplacian(Cmu_/sigmaEps_ * Tsf * Rupperf, epsilon_, "laplacian(epsilonR)");
    }

    epsEqn().relax();

    if(directImpose_ && lesCellLables_.size())
    {
        Info << "Imposing epsilonAvg on epsilonR" << endl;
        epsEqn().setValues(lesCellLables_, lesCellValuesEpsilon_);
    }

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
                - fvm::laplacian(Cmu_/sigmaK_ * Tsf * Rdiagf, k_, "laplacian(kR)")
                - fvm::laplacian(nu(), k_, "laplacian(kR)")
                ==
                G
                - fvm::Sp(epsilon_/k_, k_)
                + 0.5*tr(QR_)
            );

        if(crossTurbDiffusion_)
        {
            kEqn() -= fvc::laplacian(Cmu_/sigmaK_ * Tsf * Rupperf, k_, "laplacian(kR)");
        }
        
        kEqn().relax();

        if(directImpose_ && lesCellLables_.size())
        {
            Info << "Imposing kLes on kR" << endl;
            kEqn().setValues(lesCellLables_, lesCellValuesK_);
        }

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
            - fvm::laplacian(Cmu_/sigmaK_ * Tsf * Rdiagf, R_, "laplacian(RR)")
            - fvm::laplacian(nu(), R_, "laplacian(RR)")
            + fvm::Sp(epsilon_/k_, R_)
            ==                                        
            P                                         // production tensor
            + k_ * f_
            + QR_
        );

    if(crossTurbDiffusion_)
    {
        REqn() -= fvc::laplacian(Cmu_/sigmaK_*Ts*Rupper, R_, "laplacian(RR)");
    }

    REqn().relax();

    if(directImpose_ && lesCellLables_.size())
    {
        Info << "Imposing RAvg on RR" << endl;
        REqn().setValues(lesCellLables_, lesCellValuesR_);
    }

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

    volScalarField L_ = L();
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
            +    Cg2_ * k_/Ts * dev(symm(fbij & fbij))
            +   (Cg3_ - Cg3s_ * sqrt(bij && bij)) * k_ * Sij
            +    Cg4_ * k_ * dev( twoSymm(fbij & fSij) )
            +    Cg5_ * k_ * (twoSymm(fbij & Wij )) ;
    }
    
    
    tmp<fvSymmTensorMatrix> fEqn
        (
            fvm::laplacian(f_)
            ==
            fvm::Sp(1.0/sqr(L_), f_)
            -
            (
                exSrc/k_ + dev(R_)/(k_*Ts)
            ) / sqr(L_)
            );

    if(ellipticOperatorCorrection_)
    {
        volVectorField gradL = fvc::grad(L_);

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
        fEqn() -= fvm::Sp((16.0*beta_*magSqr(gradL))/sqr(L_), f_) - 8.0*beta_* LdF/L_;
    }
    
    fEqn().relax();
    fEqn().boundaryManipulate(f_.boundaryField());
    solve(fEqn);

    if(debug) 
    {
        Info << "  min(tr(f)) before: " << gMin(tr(f_)()) << endl;
    }
    
    if(zeroTraceF_)
    {
        f_ = dev(f_);
    }

    if(debug) 
    {
        Info << "  min(tr(f)) after zero out trace: " << gMin(tr(f_)()) << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
