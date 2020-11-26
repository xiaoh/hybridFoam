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

    ellipticOperatorCorrection_
    (
        coeffDict_.lookupOrAddDefault<Switch>
        (
            "ellipticOperatorCorrection",
            false
        )
    ),
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
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Durbin::read()
{
    if (GenElliptic::read())
    {
        solveK_.readIfPresent(word("solveK"), coeffDict());
        ellipticOperatorCorrection_.readIfPresent
        (
            word("ellipticOperatorCorrection"), coeffDict()
        );
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

    volScalarField T_("T", T());

    #include "epsilonWallI.H" // set patch internal eps values

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
        + fvm::div(phi_, epsilon_)
        - fvm::Sp(fvc::div(phi_), epsilon_)
        - fvm::laplacian(Cmu_/sigmaEps_*T_*R_ + nu()*I, epsilon_)
      ==
        C1_*G/T_*(1.0 + 0.1*G/epsilon_)
        - fvm::Sp(C2_/T_, epsilon_)
    );

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
                - fvm::laplacian(Cmu_/sigmaK_*T_*R_ + nu()*I, k_)
                ==
                G
                - fvm::Sp(epsilon_/k_, k_)
            );

        kEqn().relax();
        solve(kEqn);
    }
    else
    {
        k_ = 0.5*tr(R_);
    }
    bound(k_, k0_);

    // Reynolds stress equation
    #include "fWallI.H" // set patch internal f values

    tmp<fvSymmTensorMatrix> REqn
        (
            fvm::ddt(R_)
            + fvm::div(phi_, R_)
            - fvm::Sp(fvc::div(phi_), R_)
            - fvm::laplacian(Cmu_*T_*R_ + nu()*I, R_)
            + fvm::Sp(epsilon_/k_, R_)
            ==
            P                                         // production tensor
            + k_*f_
        );

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

    volScalarField L_ = L();
    T_ = T(); // re-compute time scale

    volSymmTensorField exSrc = -Clrr1_*dev(R_)/T_ - Clrr2_*dev(P);

    if(SSG_)
    {
        volSymmTensorField bij("bij", 0.5*dev(R_/k_));
        volTensorField gradU = fvc::grad(U_);
        volSymmTensorField Sij("Sij", dev(symm(gradU)));
        volTensorField Wij = skew(gradU);

        exSrc =
          - (Cg1_*k_/T_ + Cg1s_*G)*bij
          + Cg2_*k_/T_*dev(symm(bij & bij))
          + (Cg3_ - Cg3s_*sqrt(bij && bij))*k_*Sij
          + Cg4_*k_*dev(twoSymm(symm2full(bij) & Sij))
          + Cg5_*k_*(twoSymm(bij & Wij)) ;
    }


    tmp<fvSymmTensorMatrix> fEqn
    (
        fvm::laplacian(f_)
      ==
        fvm::Sp(1.0/sqr(L_), f_)
      - (
            exSrc/k_ + dev(R_)/(k_*T_)
        )/sqr(L_)
    );

    if(ellipticOperatorCorrection_)
    {
        volVectorField gradL = fvc::grad(L_);

        volSymmTensorField LdF
        (
            IOobject
            (
                "LnF",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedSymmTensor
            (
                "Ldf::zero",
                dimensionSet(0, -1, -1, 0, 0),
                symmTensor::zero
            )
        );

        for (label compi=0; compi<6; compi++)
        {
            volScalarField fcomp = f_.component(compi);
            LdF.replace(compi, gradL & fvc::grad(fcomp));
        }
        // Append the correction terms. (note == above is equivalent to -)
        fEqn() -=
            fvm::Sp((16.0*beta_*magSqr(gradL))/sqr(L_), f_)
          - 8.0*beta_* LdF/L_;
    }

    fEqn().relax();
    fEqn().boundaryManipulate(f_.boundaryField());
    solve(fEqn);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
