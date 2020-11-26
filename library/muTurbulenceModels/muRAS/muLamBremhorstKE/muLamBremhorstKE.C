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

#include "muLamBremhorstKE.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(muLamBremhorstKE, 0);
addToRunTimeSelectionTable(muRASModel, muLamBremhorstKE, dictionary);


// Diagnosis of difference between LES/RANS K and Epsilon quantities)
void muLamBremhorstKE::diagnosisDiffKE(string msg) const
{
    const volScalarField& lesFlag( mesh_.lookupObject<const volScalarField>("lesFlagR") );
    volScalarField kLES = scalar(0.5) * tr(RAvg_); // Update kLES (redundunt)
    volScalarField epsilonLES = epsilonAvg_;
    epsilonLES = max(epsilonLES, epsilonSmall_);

  // Set the LES/RAS difference in turbulence quantities.
    volScalarField diffKField = Foam::mag(kLES - k_);
    
    volScalarField diffEField = Foam::mag(epsilonLES - epsilon_);
    
    dimensionedScalar diffK = diffKField.weightedAverage(mesh_.V() * lesFlag);
    dimensionedScalar diffEpsilon = diffEField.weightedAverage(mesh_.V() * lesFlag);

    Info << "RANS-LES turbulence deviation : "
         << "diff(k) = " << diffK.value() << ", "
         << "diff(Epsilon) = " << diffEpsilon.value() << endl;

    dimensionedScalar kLESMask = kLES.weightedAverage(mesh_.V() * lesFlag);
    dimensionedScalar kRASMask = k_.weightedAverage(mesh_.V() * lesFlag);
    dimensionedScalar eLESMask = epsilonLES.weightedAverage(mesh_.V() * lesFlag);
    dimensionedScalar eRASMask = epsilon_.weightedAverage(mesh_.V() * lesFlag);
  
    Info << "KE on LES domain: "
         << "kLESAvg = " << kLESMask.value() << " "
         << "kRAS = "<< kRASMask.value() << " "
         << "EpsLES = " << eLESMask.value() << " "
         << "EpsRAS = " << eRASMask.value() << endl;
  
    volScalarField rasFlag = 1.0 - lesFlag;
    
    dimensionedScalar kLESMaskR = kLES.weightedAverage(mesh_.V() * rasFlag);
    dimensionedScalar kRASMaskR = k_.weightedAverage(mesh_.V() * rasFlag);
    dimensionedScalar eLESMaskR = epsilonLES.weightedAverage(mesh_.V() * rasFlag);
    dimensionedScalar eRASMaskR = epsilon_.weightedAverage(mesh_.V() * rasFlag);
    
    Info << "KE on RAS domain: "
         << "kLESAvg = " << kLESMaskR.value() << " "
         << "kRAS = "<< kRASMaskR.value() << " "
         << "EpsLES = " << eLESMaskR.value() << " "
         << "EpsRAS = " << eRASMaskR.value() << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

muLamBremhorstKE::muLamBremhorstKE
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel,
    const volSymmTensorField & RAvg,
    const volScalarField & epsilonAvg
)
:
    muRASModel(typeName, U, phi, lamTransportModel,
             RAvg, epsilonAvg),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            1.92
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaEps",
            coeffDict_,
            1.3
        )
    ),

    k_
    (
        IOobject
        (
            "kR",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    epsilon_
    (
        IOobject
        (
            "epsilonR",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    y_(mesh_),

    Rt_(sqr(k_)/(nu()*epsilon_)),

    fMu_
    (
        sqr(scalar(1) - exp(-0.0165*(sqrt(k_)*y_/nu())))
       *(scalar(1) + 20.5/(Rt_ + SMALL))
    ),

    nut_("nutR", Cmu_*fMu_*sqr(k_)/(epsilon_ + epsilonSmall_)),

  Qk_("QkR", k_/runTime_.deltaT()*0.0),
  Qepsilon_(epsilon_/runTime_.deltaT()*0.0),
  kLES_
  (
      IOobject
      (
          "mrkLESAvg",
          runTime_.timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      scalar(0.5)*tr(RAvg_)
  )
{
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> muLamBremhorstKE::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> muLamBremhorstKE::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> muLamBremhorstKE::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}


bool muLamBremhorstKE::read()
{
    if (muRASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void muLamBremhorstKE::correct()
{
    muRASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    if (mesh_.changing())
    {
        y_.correct();
    }

    volScalarField G("muRASModel::G", nut_*2*magSqr(symm(fvc::grad(U_))));


    // Calculate parameters and coefficients for low-Reynolds number model

    Rt_ = sqr(k_)/(nu()*epsilon_);
    volScalarField Ry = sqrt(k_)*y_/nu();

    fMu_ = sqr(scalar(1) - exp(-0.0165*Ry))*(scalar(1) + 20.5/(Rt_ + SMALL));

    volScalarField f1 = scalar(1) + pow(0.05/(fMu_ + SMALL), 3);
    volScalarField f2 = scalar(1) - exp(-sqr(Rt_));


    // Dissipation equation

    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        C1_*f1*G*epsilon_/k_
      - fvm::Sp(C2_*f2*epsilon_/k_, epsilon_)
      + Qepsilon_
    );

    epsEqn().relax();
    solve(epsEqn);
    bound(epsilon_, epsilon0_);


    // Turbulent kinetic energy equation

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G - fvm::Sp(epsilon_/k_, k_)
        + Qk_
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, k0_);


    // Re-calculate viscosity
    nut_ == Cmu_*fMu_*sqr(k_)/epsilon_;
}


//- enforce LES results on RAS fields via relax forcing
void muLamBremhorstKE::enforceFields
  (
   dimensionedScalar turbRelaxTime, scalar rampQ, word enforceMode
  )
  {
    const volScalarField& lesFlag( mesh_.lookupObject<const volScalarField>("lesFlagR") );
    kLES_ = scalar(0.5) * tr(RAvg_); // Update kLES first 

    volScalarField epsilonLES ("epsilonLES", epsilonAvg_);
    epsilonLES = max(epsilonLES, epsilonLES * scalar(0));

    scalar rampQLim = Foam::min(Foam::max(rampQ, scalar(0)), scalar(1));

    // bound target k and epsilon (from LES)
    bound(kLES_, scalar(2)*k0_);
    bound(epsilonLES, scalar(2)*epsilon0_);

    if(turbRelaxTime < runTime_.deltaT() * thresholdImpose_)
      {
        // Directly set fields
          Info << "Enforcing both on RAS directly." << endl;
          k_ *= (1.0 - lesFlag); // Clean out the data on LES zone
          k_ += (kLES_ * lesFlag); // Reset the data on LES zone
          epsilon_ *= (1.0 - lesFlag);
          epsilon_ += (epsilonLES * lesFlag);
          
          Qk_ *= scalar(0.0);
        Qepsilon_ *= scalar(0.0);
      }
    else
      {
        //- Use relax forcing to indirectly influnce K/E fields
        //  through the "correction" operation of next step
        // Compute the forcing on the LES zone only
          Info << "Enforcing all turb quantities on RAS via forcing." << endl;
          Qk_ = (kLES_ - k_) / turbRelaxTime * lesFlag * rampQLim;
          Qepsilon_ = 
              (epsilonLES - epsilon_) / turbRelaxTime 
              * lesFlag * rampQLim;
      }
  }

  void muLamBremhorstKE::diagnosis() const
  {
    diagnosisDiffKE(string(" Final "));   
  }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
