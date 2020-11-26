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

#include "ManceauHanjalic.H"
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

defineTypeNameAndDebug(ManceauHanjalic, 0);
addToRunTimeSelectionTable(RASModel, ManceauHanjalic, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ManceauHanjalic::ManceauHanjalic
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),
    GenElliptic(U, phi, lamTransportModel),

    alpha_
    (
        IOobject
        (
            "alpha",
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
bool ManceauHanjalic::read()
{
    if (RASModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


void ManceauHanjalic::correct()
{
    GenElliptic::correct();

    if (!turbulence_)
    {
        return;
    }

    volSymmTensorField P = -twoSymm(R_ & fvc::grad(U_));
    volScalarField G("RASModel::G", 0.5*mag(tr(P)));

    volScalarField T_("T", T());

    // Direction of alpha
    volVectorField n = fvc::grad(alpha_);
    n /= mag(n) + dimensionedScalar("nsmall", n.dimensions(), VSMALL);
    volSymmTensorField nn = symm(n*n);

    #include "epsilonWallI.H" // set patch internal eps values

    // weight of the homogeneous contribution of the pressure strain
    // indicating wall distance

    volScalarField fa = pow(alpha_, 2);

    volScalarField C1Modify = 1.0 + 0.03*(1.0 - fa)
      * sqrt
        (
            max(k_/max((R_ && nn), k0_), VSMALL)
        );

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::Sp(fvc::div(phi_), epsilon_)
      - fvm::laplacian(Cmu_/sigmaEps_*T_*R_ + nu()*I, epsilon_)
   // - fvm::laplacian(DepsilonEff(), epsilon_)
      ==
        C1_*G/T_*C1Modify
      - fvm::Sp(C2_/T_, epsilon_)
    );

    epsEqn().relax();
    epsEqn().boundaryManipulate(epsilon_.boundaryField());

    solve(epsEqn);
    bound(epsilon_, epsilon0_);

    // wall contribution of pressure strain (weight included). Note:
    // repeated typo in Eqs. (9) and (A.10) in (Manceau & Hanjalic 2002)
    volSymmTensorField psw =
      - 5.0*(1.0 - fa)*epsilon_/k_
      * (
            twoSymm(R_ & nn) - 0.5*(R_ && nn)*(nn + I)
        );

    // used LRR-IP as default; SSG is also option.
    volScalarField imSrc = Clrr1_*fa*epsilon_/k_;
    volSymmTensorField exSrc =
        (Clrr1_*twoThirdsI)*fa*epsilon_ - fa*Clrr2_*dev(P);

    if (SSG_) // SSG model. Disabled by default
    {
        volSymmTensorField bij("bij", 0.5*dev(R_/k_));
        volTensorField gradU = fvc::grad(U_);
        volSymmTensorField Sij("Sij", dev(symm(gradU)));
        volTensorField Wij = skew(gradU);

        imSrc = 0.5*fa*(Cg1_/T_ + Cg1s_*G/k_);
        exSrc = fa
          * (
                (Cg1_*k_/T_ + Cg1s_*G)*oneThirdI
              + Cg2_*k_/T_*dev(symm(bij & bij))
              + (Cg3_ - Cg3s_*sqrt(bij && bij))*k_*Sij
              + Cg4_*k_*dev(twoSymm(symm2full(bij) & Sij))
              + Cg5_*k_*twoSymm(bij & Wij)
            );
    }

    // Reynolds stress equation
    tmp<fvSymmTensorMatrix> REqn
        (
            fvm::ddt(R_)
          + fvm::div(phi_, R_)
          - fvm::Sp(fvc::div(phi_), R_)
          - fvm::laplacian(Cmu_*T_*R_ + nu()*I, R_) // Daly-Harlow diffusion model
      //  - fvm::laplacian(DREff(), R_)       // isotropized diffusion
          + fvm::Sp(imSrc, R_)                // pressure-rate-of-strain, implicit (RHS) source
          + fvm::Sp((1 - fa)*epsilon_/k_, R_) // wall contribution of epsilon
          ==                                  // ==
            P                                 // production tensor
          - twoThirdsI*epsilon_*fa            // homogeneous contribution of epsilon
          + exSrc                             // pressure-rate-of-strain, explicit (LHS) source
          + psw                               // pressure-rate-of-strain, wall contribution
        );

    REqn().relax();
    solve(REqn);

    R_.max
    (
        dimensionedSymmTensor
        (
            "zero",
            R_.dimensions(),
            symmTensor
            (
                k0_.value(), -GREAT, -GREAT,
                k0_.value(), -GREAT,
                k0_.value()
            )
       )
    );

    k_ = 0.5*tr(R_);

    bound(k_, k0_);

    volScalarField L_ = L();
    T_ = T(); // re-compute time scale

    tmp<fvScalarMatrix> alphaEqn
    (
        fvm::laplacian(alpha_)
     ==
        fvm::Sp(1.0/sqr(L_), alpha_)
      - scalar(1.0)/(sqr(L_))
    );

    alphaEqn().relax();
    solve(alphaEqn);

    // Finally, re-calculate turbulent viscosity according to V2F model
    // nut_ = 2.0/3.0*Cmu_*max((R_ && nn), k0_)*T_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
