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

#include "muKOmega.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(muKOmega, 0);
addToRunTimeSelectionTable(muRASModel, muKOmega, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Diagnosis of difference between LES/RANS K and Epsilon quantities)
void muKOmega::diagnosisDiffKE(string msg) const
{
    if(turbConsistencyDiagnosis_)
    {
        const volScalarField& lesFlag( mesh_.lookupObject<const volScalarField>("lesFlagR") );
        dimensionedScalar vsmall("vsmall_lau", dimVolume, SMALL);
        volScalarField omegaLES ("omegaLES", epsilonAvg_/(Cmu_ * k_));
        omegaLES = max(omegaLES, omega0_);
        
        // Set the LES/RAS difference in turbulence quantities.
        volScalarField diffKField = Foam::mag(kLES_ - k_);
        
        volScalarField diffEField = Foam::mag(omegaLES - omega_);
        
        dimensionedScalar diffK = diffKField.weightedAverage(mesh_.V() * lesFlag + vsmall);
        dimensionedScalar diffOmega = diffEField.weightedAverage(mesh_.V() * lesFlag + vsmall);
        
        Info << "RANS-LES turbulence deviation : "
             << "diff(k) = " << diffK.value() << ", "
             << "diff(Omega) = " << diffOmega.value() << endl;
        
        dimensionedScalar kLESMask = kLES_.weightedAverage(mesh_.V() * lesFlag + vsmall);
        dimensionedScalar kRASMask = k_.weightedAverage(mesh_.V() * lesFlag + vsmall);
        dimensionedScalar eLESMask = omegaLES.weightedAverage(mesh_.V() * lesFlag + vsmall);
        dimensionedScalar eRASMask = omega_.weightedAverage(mesh_.V() * lesFlag + vsmall);
        
        Info << "KE on LES domain: "
             << "kLESAvg = " << kLESMask.value() << " "
             << "kRAS = "<< kRASMask.value() << " "
             << "OmgLES = " << eLESMask.value() << " "
             << "OmeRAS = " << eRASMask.value() << endl;
        
        volScalarField rasFlag = 1.0 - lesFlag;
        
        dimensionedScalar kLESMaskR = kLES_.weightedAverage(mesh_.V() * rasFlag + vsmall);
        dimensionedScalar kRASMaskR = k_.weightedAverage(mesh_.V() * rasFlag + vsmall);
        dimensionedScalar eLESMaskR = omegaLES.weightedAverage(mesh_.V() * rasFlag + vsmall);
        dimensionedScalar eRASMaskR = omega_.weightedAverage(mesh_.V() * rasFlag + vsmall);
        
        Info << "KE on RAS domain: "
             << "kLESAvg = " << kLESMaskR.value() << " "
             << "kRAS = "<< kRASMaskR.value() << " "
             << "EpsLES = " << eLESMaskR.value() << " "
             << "EpsRAS = " << eRASMaskR.value() << endl;
    }
}


void muKOmega::checkEnforceMode(word enforceMode)
{

  if (
      enforceMode != "kOnly" 
      && enforceMode != "both"
      && enforceMode != "none"
      )
    {
      FatalErrorIn("muKOmega::checkEnforceMode()")
        << "enforceMode should be one of the following: both, kOnly, none!\n"
        << "Instead you have: " 
        << enforceMode
        << abort(FatalError);
    }
}


muKOmega::muKOmega
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
            "betaStar",
            coeffDict_,
            0.09
        )
    ),
    beta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta",
            coeffDict_,
            0.072
        )
    ),
    alpha_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alpha",
            coeffDict_,
            0.52
        )
    ),
    alphaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK",
            coeffDict_,
            0.5
        )
    ),
    alphaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega",
            coeffDict_,
            0.5
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
    omega_
    (
        IOobject
        (
            "omegaR",
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
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
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

  Qomega_
    (
        IOobject
        (
            "QOmega",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qOmg", dimensionSet(0, 0, -2, 0, 0), 0)
    ),

  kLES_
  (
      IOobject
      (
          "mrkLESAvg",
          runTime_.timeName(),
          mesh_,
          // IOobject::NO_READ,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
      ),
      scalar(0.5)*tr(RAvg_)
  ),

    directImpose_(lookupOrDefault<Switch>("directImpose", false)),
    lesCellLables_(),
    lesCellValuesK_(),
    lesCellValuesOmega_(),
    imposeTurbEvery_(1),
    omegaMax_(GREAT)

{
    nut_ = k_/(omega_ + omegaSmall_);
    nut_.correctBoundaryConditions();

    IOdictionary relaxParameters
        (
            IOobject
            (
                "relaxParameters",
                runTime_.constant(),
                "../../constant",
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

    dictionary couplingDict(relaxParameters.subDictPtr("couplingOptions"));

    imposeTurbEvery_ = 
        couplingDict.lookupOrDefault<label>("mapL2REvery", 1, true);

    omegaMax_ = couplingDict.lookupOrDefault<scalar>("omegaMax", GREAT, true);
    
    Info << "(If enabled) directly imposing turb quantities every: " << imposeTurbEvery_ << endl;
        
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> muKOmega::R() const
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


tmp<volSymmTensorField> muKOmega::devReff() const
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


tmp<fvVectorMatrix> muKOmega::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}


bool muKOmega::read()
{
    if (muRASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        beta_.readIfPresent(coeffDict());
        alphaK_.readIfPresent(coeffDict());
        alphaOmega_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void muKOmega::correct()
{
    muRASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    volScalarField G("muRASModel::G", nut_*2*magSqr(symm(fvc::grad(U_))));

    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::Sp(fvc::div(phi_), omega_)
      - fvm::laplacian(DomegaEff(), omega_)
     ==
        alpha_*G*omega_/k_
        - fvm::Sp(beta_*omega_, omega_)
        + fvm::SuSp(Qomega_/omega_, omega_)
    );

    omegaEqn().relax();

    if(directImpose_ && lesCellLables_.size())
    {
        Info << "Imposing omegaLES on omegaR" << endl;
        omegaEqn().setValues(lesCellLables_, lesCellValuesOmega_);
    }
    
    omegaEqn().boundaryManipulate(omega_.boundaryField());

    solve(omegaEqn);
    bound(omega_, omega0_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::Sp(fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_)
      ==
        G
        - fvm::Sp(Cmu_*omega_, k_)
        + fvm::SuSp(Qk_/k_, k_)
    );

    kEqn().relax();

    if(directImpose_ && lesCellLables_.size())
    {
        Info << "Imposing kLES on kR" << endl;
        kEqn().setValues(lesCellLables_, lesCellValuesK_);
    }

    solve(kEqn);
    bound(k_, k0_);


    // Re-calculate viscosity
    nut_ = k_/(omega_ + omegaSmall_);
    nut_.correctBoundaryConditions();
}


//- enforce LES results on RAS fields via relax forcing
  void muKOmega::enforceFields
  (
   dimensionedScalar turbRelaxTime, scalar rampQ, word enforceMode
  )
{
    const volScalarField& lesFlag( mesh_.lookupObject<const volScalarField>("lesFlagR") );
    kLES_ = scalar(0.5) * tr(RAvg_); // Update kLES first 
    volScalarField omegaLES ("omegaLES", epsilonAvg_/(Cmu_ * k_));
    label timeIndex = runTime_.timeIndex();

    scalar rampQLim = Foam::min(Foam::max(rampQ, scalar(0)), scalar(1));

    checkEnforceMode(enforceMode);

    // bound target k and omega (from LES)
    kLES_ = max( min(kLES_, k_ * boundFactor_), k_ / boundFactor_);
    omegaLES = max( min(omegaLES, omega_ * boundFactor_), omega_ / boundFactor_);

    if(turbRelaxTime < runTime_.deltaT() * thresholdImpose_)
      {
          // prepare lists for LES cell labels and values 
          // (for setting values in corresponding matrix later)
          if( timeIndex % imposeTurbEvery_ == 0)
          {
              directImpose_ = true;
              Info << "Update imposed values." << endl;
              label N = sum(pos(lesFlag.internalField()-0.5));
              
              lesCellLables_.clear(); 
              lesCellValuesK_.clear();
              lesCellValuesOmega_.clear();
              
              lesCellLables_.setSize(N); 
              lesCellValuesK_.setSize(N);
              lesCellValuesOmega_.setSize(N);
              
              label i=0;
              forAll(lesFlag, celli) // scan lesFlag Field
              {
                  if(lesFlag[celli] > 0.5)
                  {
                      lesCellLables_[i] = celli;
                      lesCellValuesK_[i] = kLES_[celli];
                      lesCellValuesOmega_[i] = omegaLES[celli];
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
          
          Qk_ *= scalar(0.0);
          Qomega_ *= scalar(0.0);
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
            Qomega_ *= scalar(0.0);
          }
        else if (enforceMode == "both") 
          {
            Info << "Enforcing all turb quantities on RAS via forcing." << endl;
            Qk_ = (kLES_ - k_) / turbRelaxTime * lesFlag * rampQLim;
            Qomega_ = 
              (omegaLES - omega_) / turbRelaxTime 
              * lesFlag * rampQLim;
          }
        else if (enforceMode == "none")  
        {
            Info << "Not enforcing any turb quantities on RAS." << endl;
            Qk_ *= scalar(0.0);
            Qomega_ *= scalar(0.0);
        }
        else
        {
            Info << "muKOmega::enforceFields:: something went wrong with the enforceMode.\n"
                 << "Should quit now." << endl;
        }
      }
  }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
