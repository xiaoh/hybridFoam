/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

#include "muV2F.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(muV2F, 0);
addToRunTimeSelectionTable(muRASModel, muV2F, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> muV2F::T() const
{
    volScalarField yStar_=pow(CmuKE,0.25)*sqrt(k_)*yw_/nu();
    return max
        (
            k_/(epsilon_ + epsilonSmall_),
            pos(yStarLim - yStar_) * 6.0 * sqrt(nu()/(epsilon_ + epsilonSmall_))
        );
}

tmp<volScalarField> muV2F::L() const
{
    volScalarField yStar_=pow(CmuKE,0.25)*sqrt(k_)*yw_/nu();
    return
        CL*max
        (
            pow(k_, 1.5) / (epsilon_ + epsilonSmall_),
            pos(yStarLim-yStar_) * CEta * pow(pow(nu(),3.0)/(epsilon_ + epsilonSmall_),0.25)
        );
}


// Diagnosis of difference between LES/RANS K and Epsilon quantities)
void muV2F::diagnosisDiffKE(string msg) const
{
    const volScalarField& lesFlag( mesh_.lookupObject<const volScalarField>("lesFlagR") );
    volScalarField kLES = scalar(0.5) * tr(RAvg_); // Update kLES (redundunt)
    const volScalarField& epsilonLES = epsilonAvg_; //data member of muTurbulenceModel
    dimensionedScalar vsmall("vsmall_v2f", dimVolume, SMALL);

    // Set the LES/RAS difference in turbulence quantities.
    volScalarField diffKField = Foam::mag(kLES - k_);
    
    volScalarField diffEField = Foam::mag(epsilonLES - epsilon_);
    
    dimensionedScalar diffK = diffKField.weightedAverage(mesh_.V() * lesFlag + vsmall);
    dimensionedScalar diffEpsilon = diffEField.weightedAverage(mesh_.V() * lesFlag + vsmall);
    
    Info << "RANS-LES turbulence deviation : "
         << "diff(k) = " << diffK.value() << ", "
         << "diff(Epsilon) = " << diffEpsilon.value() << endl;
    
    dimensionedScalar kLESMask = kLES.weightedAverage(mesh_.V() * lesFlag + vsmall);
    dimensionedScalar kRASMask = k_.weightedAverage(mesh_.V() * lesFlag + vsmall);
    dimensionedScalar eLESMask = epsilonLES.weightedAverage(mesh_.V() * lesFlag + vsmall);
    dimensionedScalar eRASMask = epsilon_.weightedAverage(mesh_.V() * lesFlag + vsmall);
    
    Info << "KE on LES domain: "
         << "kLESAvg = " << kLESMask.value() << " "
         << "kRAS = "<< kRASMask.value() << " "
         << "EpsLES = " << eLESMask.value() << " "
         << "EpsRAS = " << eRASMask.value() << endl;
    
    volScalarField rasFlag = 1.0 - lesFlag;
    
    dimensionedScalar kLESMaskR = kLES.weightedAverage(mesh_.V() * rasFlag + vsmall);
    dimensionedScalar kRASMaskR = k_.weightedAverage(mesh_.V() * rasFlag + vsmall);
    dimensionedScalar eLESMaskR = epsilonLES.weightedAverage(mesh_.V() * rasFlag + vsmall);
    dimensionedScalar eRASMaskR = epsilon_.weightedAverage(mesh_.V() * rasFlag + vsmall);
    
    Info << "KE on RAS domain: "
         << "kLESAvg = " << kLESMaskR.value() << " "
         << "kRAS = "<< kRASMaskR.value() << " "
         << "EpsLES = " << eLESMaskR.value() << " "
         << "EpsRAS = " << eRASMaskR.value() << endl;
}


void muV2F::checkEnforceMode(word enforceMode)
{

  if (
      enforceMode != "kOnly" 
      && enforceMode != "both"
      && enforceMode != "none"
      )
    {
      FatalErrorIn("muV2F::checkEnforceMode()")
        << "enforceMode should be one of the following: both, kOnly, none!\n"
        << "Instead you have: " 
        << enforceMode
        << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
muV2F::muV2F
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

    Cmu
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.22
        )
    ),
    
    
    CmuKE
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CmuKE",
            coeffDict_,
            0.09
        )
    ),

    Ceps10
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps10",
            coeffDict_,
            1.4
        )
    ),
    
    Ceps11
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps11",
            coeffDict_,
            0.05
        )
    ),
    
    Ceps2
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps2",
            coeffDict_,
            1.9
        )
    ),
    
    C1
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.4
        )
    ),
    C2
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            0.3
        )
    ),
    CL
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CL",
            coeffDict_,
            0.23
        )
    ),
    CEta
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CEta",
            coeffDict_,
            70.0
        )
    ),
    
    oneOnSigmaK
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "oneOnSigmaK",
            coeffDict_,
            1.0
        )
    ),
    
    oneOnSigmaEps
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "oneOnSigmaEps",
            coeffDict_,
            0.77
        )
    ),
    
    yStarLim
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "yStarLim",
            coeffDict_,
            5.0
        )
    ),
    

    f0_("f0small", dimless/dimTime, SMALL),

    yw_(mesh_),

    k_
    (
        IOobject
        (
            "k",
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
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    v2_
    (
        IOobject
        (
            "v2",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
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
    ),

    // Calculate viscosity (with Davidson correction - 2003)
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
        min
        (
            CmuKE*sqr(k_)/(epsilon_ + epsilonSmall_),
            Cmu*v2_*T()
        )
    ),

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

tmp<volSymmTensorField> muV2F::R() const
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
            ((2.0/3.0)*I)*k_ - nut_*2*symm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}

tmp<volSymmTensorField> muV2F::devReff() const
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

tmp<fvVectorMatrix> muV2F::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}


bool muV2F::read()
{
    if (muRASModel::read())
    {
        Cmu.readIfPresent(coeffDict_);
        CmuKE.readIfPresent(coeffDict_);
        Ceps10.readIfPresent(coeffDict_);
        Ceps11.readIfPresent(coeffDict_);
        Ceps2.readIfPresent(coeffDict_);
        C1.readIfPresent(coeffDict_);
        C2.readIfPresent(coeffDict_);
        CL.readIfPresent(coeffDict_);
        CEta.readIfPresent(coeffDict_);
        oneOnSigmaK.readIfPresent(coeffDict_);
        oneOnSigmaEps.readIfPresent(coeffDict_);
        yStarLim.readIfPresent(coeffDict_);
        
        return true;
    }
    else
    {
        return false;
    }
}


void muV2F::correct()
{
    transportModel_.correct();

    if (!turbulence_)
    {
        return;
    }

    muRASModel::correct();

    if (mesh_.moving())
    {
        yw_.correct();
    }

    volScalarField S2 = 2*magSqr(symm(fvc::grad(U_)));

    volScalarField G = nut_*S2;

    volScalarField T_ = T();
    volScalarField Ceps1 = Ceps10*(scalar(1.0)+Ceps11*min(sqrt(k_/v2_),scalar(10.0)));

    // Dissipation rate equation

    #   include "muEpsilonV2FWallI.H"
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        Ceps1*G/T_
        - fvm::Sp(Ceps2/T_, epsilon_)
        + Qepsilon_
    );

    epsEqn().relax();
#   include "muWallDissipationV2FI.H"
    solve(epsEqn);
    bound(epsilon_, epsilon0_);


    // Turbulent kinetic energy equation

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G - fvm::Sp(1.0/T_, k_)
//      G - fvm::Sp(epsilon_/k_, k_)
        + Qk_
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, k0_);

    // f equation
    volScalarField L_ = L();

    tmp<fvScalarMatrix> fEqn
    (
//     - fvm::laplacian(scalar(1.0),f_)
      - fvm::laplacian(f_)
     ==
      - fvm::Sp(1.0/sqr(L_),f_)
      - ((C1-scalar(6.0))*v2_/k_ - 2.0/3.0*(C1-scalar(1.0)))/(sqr(L_)*T_)
      + C2*G/(k_*sqr(L_))
      /*      - fvm::Sp(1.0/sqr(L_),f_)
          - ((C1-scalar(1.0))*(v2_/k_-scalar(2.0/3.0))/T_
          - C2*G/k_
          - 5.0*epsilon_*v2_/sqr(k_))/sqr(L_)  */
    );

    fEqn().relax();
    solve(fEqn);
    bound(f_, f0_);

    // v2 equation

    tmp<fvScalarMatrix> v2Eqn
    (
        fvm::ddt(v2_)
      + fvm::div(phi_, v2_)
      - fvm::laplacian(DkEff(), v2_)
     ==
        //    k_ * f_ 
        // Davidson correction - 2003
        min
        (
            k_ * f_, 
            -((C1-scalar(6.0))*v2_ - 2.0/3.0*k_*(C1-scalar(1.0)))/T_ + C2*G
        )
        - fvm::Sp(6.0*epsilon_/k_, v2_)
    );
    
    v2Eqn().relax();
    solve(v2Eqn);
    bound(v2_, k0_);

    // Re-calculate viscosity (with Davidson correction - 2003)
    nut_ = min
        (
            CmuKE*sqr(k_)/(epsilon_ + epsilonSmall_),
            Cmu*v2_*T()
        );
}


//- enforce LES results on RAS fields via relax forcing
  void muV2F::enforceFields
  (
    dimensionedScalar turbRelaxTime, scalar rampQ, word enforceMode
  )
  {
    const volScalarField& lesFlag( mesh_.lookupObject<const volScalarField>("lesFlag") );
    kLES_ = scalar(0.5) * tr(RAvg_); // Update kLES first 
    const volScalarField& epsilonLES = epsilonAvg_;

    scalar rampQLim = Foam::min(Foam::max(rampQ, scalar(0)), scalar(1));

    checkEnforceMode(enforceMode);

    // bound target k and epsilon (from LES)
    bound(kLES_, scalar(2)*k0_);

    if(turbRelaxTime < runTime_.deltaT() * thresholdImpose_)
      {
        // Directly set fields
        if (enforceMode == "kOnly")
          {
            Info << "Enforcing k on RAS directly." << endl;
            k_ *= (1.0 - lesFlag); // Clean out the data on LES zone
            k_ += (kLES_ * lesFlag); // Reset the data on LES zone
          }
        else if (enforceMode == "both")
          {
            Info << "Enforcing both on RAS directly." << endl;
            k_ *= (1.0 - lesFlag); // Clean out the data on LES zone
            k_ += (kLES_ * lesFlag); // Reset the data on LES zone
            epsilon_ *= (1.0 - lesFlag);
            epsilon_ += (epsilonLES * lesFlag);
          }
        else if (enforceMode == "none")
          {
            Info << "Not directly enforcing KE on RAS." << endl;
          }
        else
          {
            Info << "muV2F::enforceFields:: something went wrong with the enforceMode.\n"
                 << "Should quit now." << endl;
          }

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
              (epsilonLES - epsilon_) / turbRelaxTime 
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
            Info << "muV2F::enforceFields:: something went wrong with the enforceMode.\n"
                 << "Should quit now." << endl;
          }
      }
  }

  void muV2F::diagnosis() const
  {
    diagnosisDiffKE(string(" Final "));   
  }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
