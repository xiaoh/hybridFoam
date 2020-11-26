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

#include "muLaunderSharmaKE.H"
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

defineTypeNameAndDebug(muLaunderSharmaKE, 0);
addToRunTimeSelectionTable(muRASModel, muLaunderSharmaKE, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> muLaunderSharmaKE::fMu() const
{
    return exp(-3.4/sqr(scalar(1) + sqr(k_)/(nu()*epsilonTilda_)/50.0));
}


tmp<volScalarField> muLaunderSharmaKE::f2() const
{
    return
        scalar(1)
      - 0.3*exp(-min(sqr(sqr(k_)/(nu()*epsilonTilda_)), scalar(50.0)));
}

// Diagnosis of difference between LES/RANS K and Epsilon quantities)
void muLaunderSharmaKE::diagnosisDiffKE(string msg) const
{
    if(turbConsistencyDiagnosis_)
    {
        const volScalarField& lesFlag( mesh_.lookupObject<const volScalarField>("lesFlagR") );
        volScalarField kLES = scalar(0.5) * tr(RAvg_); // Update kLES (redundunt)
        volScalarField DLES = 2.0*nu()*magSqr(fvc::grad(sqrt(mag(kLES))));
        volScalarField epsilonTildaLES = epsilonAvg_ - DLES;
        dimensionedScalar vsmall("vsmall_lau", dimVolume, SMALL);
        epsilonTildaLES = max(epsilonTildaLES, epsilonTildaLES * scalar(0));
        
        // Set the LES/RAS difference in turbulence quantities.
        volScalarField diffKField = Foam::mag(kLES - k_);
        
        volScalarField diffEField = Foam::mag(epsilonTildaLES - epsilonTilda_);
        
        dimensionedScalar diffK = diffKField.weightedAverage(mesh_.V() * lesFlag + vsmall);
        dimensionedScalar diffEpsilon = diffEField.weightedAverage(mesh_.V() * lesFlag + vsmall);
        
        Info << "RANS-LES turbulence deviation : "
             << "diff(k) = " << diffK.value() << ", "
             << "diff(Epsilon) = " << diffEpsilon.value() << endl;
        
        dimensionedScalar kLESMask = kLES.weightedAverage(mesh_.V() * lesFlag + vsmall);
        dimensionedScalar kRASMask = k_.weightedAverage(mesh_.V() * lesFlag + vsmall);
        dimensionedScalar eLESMask = epsilonTildaLES.weightedAverage(mesh_.V() * lesFlag + vsmall);
        dimensionedScalar eRASMask = epsilonTilda_.weightedAverage(mesh_.V() * lesFlag + vsmall);
        
        Info << "KE on LES domain: "
             << "kLESAvg = " << kLESMask.value() << " "
             << "kRAS = "<< kRASMask.value() << " "
             << "EpsLES = " << eLESMask.value() << " "
             << "EpsRAS = " << eRASMask.value() << endl;
        
        volScalarField rasFlag = 1.0 - lesFlag;
        
        dimensionedScalar kLESMaskR = kLES.weightedAverage(mesh_.V() * rasFlag + vsmall);
        dimensionedScalar kRASMaskR = k_.weightedAverage(mesh_.V() * rasFlag + vsmall);
        dimensionedScalar eLESMaskR = epsilonTildaLES.weightedAverage(mesh_.V() * rasFlag + vsmall);
        dimensionedScalar eRASMaskR = epsilonTilda_.weightedAverage(mesh_.V() * rasFlag + vsmall);
        
        Info << "KE on RAS domain: "
             << "kLESAvg = " << kLESMaskR.value() << " "
             << "kRAS = "<< kRASMaskR.value() << " "
             << "EpsLES = " << eLESMaskR.value() << " "
             << "EpsRAS = " << eRASMaskR.value() << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

muLaunderSharmaKE::muLaunderSharmaKE
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
            "sigmaEps",
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

    epsilonTilda_
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

    nut_
    (
        IOobject
        (
            "nutR",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nutR", mesh_)
    ),

  Qk_
    (
        IOobject
        (
            "QkR",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qk", dimVelocity*dimVelocity/dimTime, 0)
    ),

  Qepsilon_
    (
        IOobject
        (
            "QEpsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qEps", dimensionSet(0, 2, -4, 0, 0), 0)
    ),

  kLES_
  (
      IOobject
      (
          "mrkLESAvg",
          runTime_.timeName(),
          mesh_,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE
      ),
      scalar(0.5)*tr(RAvg_)
  )

{
    nut_ = Cmu_*fMu()*sqr(k_)/(epsilonTilda_ + epsilonSmall_);
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> muLaunderSharmaKE::R() const
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


tmp<volSymmTensorField> muLaunderSharmaKE::devReff() const
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


tmp<fvVectorMatrix> muLaunderSharmaKE::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}


bool muLaunderSharmaKE::read()
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


void muLaunderSharmaKE::correct()
{
    muRASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    volScalarField S2 = 2*magSqr(symm(fvc::grad(U_)));

    volScalarField G = nut_*S2;

    volScalarField E = 2.0*nu()*nut_*fvc::magSqrGradGrad(U_);
    
    volScalarField D = 2.0*nu()*magSqr(fvc::grad(sqrt(k_)));

    // Dissipation rate equation

    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilonTilda_)
      + fvm::div(phi_, epsilonTilda_)
      - fvm::laplacian(DepsilonEff(), epsilonTilda_)
     ==
        C1_*G*epsilonTilda_/k_
        - fvm::Sp(C2_*f2()*epsilonTilda_/k_, epsilonTilda_)
        + E
        + Qepsilon_
    );

    epsEqn().relax();
    solve(epsEqn);
    bound(epsilonTilda_, epsilon0_);

    // Turbulent kinetic energy equation

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G - fvm::Sp((epsilonTilda_ + D)/k_, k_)
        + Qk_
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, k0_);

    // Re-calculate viscosity
    nut_ = Cmu_*fMu()*sqr(k_)/epsilonTilda_;
    nut_.correctBoundaryConditions();
}


//- enforce LES results on RAS fields via relax forcing
  void muLaunderSharmaKE::enforceFields
  (
   dimensionedScalar turbRelaxTime, scalar rampQ, word enforceMode
  )
  {
    const volScalarField& lesFlag( mesh_.lookupObject<const volScalarField>("lesFlagR") );
    kLES_ = scalar(0.5) * tr(RAvg_); // Update kLES first 
    volScalarField DLES = 2.0*nu()*magSqr(fvc::grad(sqrt(mag(kLES_))));
    volScalarField epsilonTildaLES ("epsilonTildaLES", epsilonAvg_ - DLES);

    scalar rampQLim = Foam::min(Foam::max(rampQ, scalar(0)), scalar(1));

    // bound target k and epsilon (from LES)
    kLES_ = max( min(kLES_, k_ * boundFactor_), k_ / boundFactor_);
    epsilonTildaLES = max( min(epsilonTildaLES, epsilonTilda_ * boundFactor_), epsilonTilda_ / boundFactor_);

    if(turbRelaxTime < runTime_.deltaT() * thresholdImpose_)
      {
        // Directly set fields
          Info << "Enforcing both on RAS directly." << endl;
          k_ *= (1.0 - lesFlag); // Clean out the data on LES zone
          k_ += (kLES_ * lesFlag); // Reset the data on LES zone
          epsilonTilda_ *= (1.0 - lesFlag);
          epsilonTilda_ += (epsilonTildaLES * lesFlag);
          
          Qk_ *= scalar(0.0);
          Qepsilon_ *= scalar(0.0);
      }
    else
      {
        //- Use relax forcing to indirectly influnce K/E fields
        //  through the "correction" operation of next step
        // Compute the forcing on the LES zone only
          if (enforceMode == "kOnly")
          { 
              Info << "Only enforcing LES k on RAS." << endl;
              Qk_ = (kLES_ - k_) / turbRelaxTime * lesFlag * rampQLim;
              Qepsilon_ *= scalar(0.0);
          }
          else if (enforceMode == "both") 
          {
              Info << "Enforcing all turb quantities on RAS via forcing." << endl;
              Qk_ = (kLES_ - k_) / turbRelaxTime * lesFlag * rampQLim;
              Qepsilon_ = 
                  (epsilonTildaLES - epsilonTilda_) / turbRelaxTime 
                  * lesFlag * rampQLim;
          }
          else if (enforceMode == "none")  
          {
              Info << "Not enforcing any turb quantities on RAS." << endl;
              Qk_ *= scalar(0.0);
              Qepsilon_ *= scalar(0.0);
          }
          else
          {
              Info << "muLaunderSharmaKE::enforceFields:: something went wrong with the enforceMode.\n"
                   << "Should quit now." << endl;
          }
      }
  }

  void muLaunderSharmaKE::diagnosis() const
  {
    diagnosisDiffKE(string(" Final "));   
  }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
