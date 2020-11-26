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

  Class
  relaxForcing

  Description
  A class for the relaxation forcing in LES/RANS coupling

  SourceFiles
  relaxForcing.H
  relaxForcing.C


  \*---------------------------------------------------------------------------*/

#ifndef relaxForcing_H
#define relaxForcing_H

#include "Switch.H"
#include "muTurbulenceModel.H"
#include "IOdictionary.H"
#include "muRASModel.H"
#include "muLESModel.H"
#include "meshTalkChannel.H"
#include "fvCFD.H"
#include "error.H"
#include "objectRegistry.H"
#include "resolutionModel.H"
#include "staticModel.H"
#include "wallVicinity.H"
#include "Random.H"
#include "labelList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  /*---------------------------------------------------------------------------* \
    Class relaxForcing Declaration
  \*---------------------------------------------------------------------------*/


  class relaxForcing
  {
            
    // Private data
  private:

      //- References to velocity and pressure fields
      //  and turbulence models
      const volVectorField& UL_;
      const volVectorField& UR_;

      incompressible::muLESModel& sgsModelRef_;
      incompressible::muRASModel& ransModelRef_;

      // Consistent velocity field
      volVectorField& UAvg_;
      // Consistent Reynolds stress <u_i u_j>^LES
      volSymmTensorField& RAvg_;
      // Field for total kinetic energy
      volSymmTensorField& RResolvedAvg_;
      // Field for total sgs energy
      volSymmTensorField& RSgsAvg_;
      // Consistent dissipation field
      volScalarField& epsilonAvg_;

      // Forcing for LES
      volVectorField& QLes_;
      // Forcing for RANS
      volVectorField& QRans_;

      // References to the database and meshes
      const Time& runTime_;
      const fvMesh& meshL_;
      const fvMesh& meshR_;
      const label nCellsGL_;
      const label nCellsGR_;
      
      // Dictionaries
      IOdictionary & relaxParameters_;
      dictionary     timeScaleDict_;
      dictionary     couplingDict_;
      
      // SGS modeling in terms of forcing on filtered momentum equation
      volVectorField QSgs_;
      volVectorField QSgsAvg_;

      // Time scale for averaging (consistent quantities)
      const dimensionedScalar averagingTime_;
      // Time scale for relaxation
      const dimensionedScalar URelaxTime_;
      const dimensionedScalar RRelaxTime_;
      const dimensionedScalar MRelaxTime_;
      const dimensionedScalar turbRelaxTime_;
      const dimensionedScalar nu0_;
      const dimensionedScalar nuRelaxFactor_;
      const dimensionedScalar nuRelaxFactorL_;
      const dimensionedScalar nuRelaxFactorR_;
      // Whether to enforce LES fields on RAS, and
      // whether to add Q to the momentum equations
      const Switch LesRansCoupling_;
      word  RConsisLevel_;
      symmTensor flexMask_;
      word  epsilonSgsBy_;
      Switch rmDivTauSgs_;
      const dimensionedScalar couplingStartTime_;
      const dimensionedScalar gradualCouplingDuration_;
      IOdictionary consistentDict_;
      scalar rampQ_;
      word enforceMode_;
      const label mapR2LEvery_;
      const label mapL2REvery_;
      
      Switch couplingStatus_;
      Switch neutralQMomentum_;
      
      const dimensionedScalar ksmall_;
      const dimensionedScalar vsmall_;
      // Averaged SGS kinetic energy (come after mesh etc)
      volScalarField kSgsAvg_;

      volScalarField kAvg_;
      volScalarField kLES_; // total instantaneous LES energy, output only
      autoPtr<volSymmTensorField> RCorrPtr_; // Correction stress for direct forcing

      autoPtr<resolutionModel> turbFlag_;
      const volScalarField& lesFlagL_;
      const volScalarField& lesFlagR_;
      // meshTalkChannel (for the L/R meshes to communicate)
      meshTalkChannel meshTalk_;

      Switch consistencyDiagnosis_;
      Switch forcingDiagnosis_;

    // Private Member Functions

    //- Update integral length scale and flag
    void evaluateResolution();

    //- Update consistent quantities
    void updateConsistent();

    //- Compute forcing terms QLES/QRans
    void computeForcing();

    void initConsistentDict();

    //check fields are presents if restart; re-initialize fields if fresh start.
    void checkRestart();

    void updateCouplingStatus();

    //- Diagnosis of difference of UL & UR
    void diagnosisDiffU();

    //- Print relaxation parameters
    void printRelaxParameters();

    void checkRCorrPtr();
      
    //- Disallow default bitwise copy construct
      // relaxForcing(const relaxForcing&);

    //- Disallow default bitwise assignment
      // void operator=(const relaxForcing&);


  public:


    // Constructors

    //- Construct from components
    relaxForcing
    (            
        const volVectorField& UL,
        const volVectorField& UR,
        Foam::incompressible::muLESModel& sgsModelRef,
        Foam::incompressible::muRASModel& ransModelRef,            
        volVectorField& UAvg,
        volSymmTensorField& RAvg,
        volSymmTensorField& RResolvedAvg,
        volSymmTensorField& RSgsAvg,
        volScalarField& epsilonAvg,
        volVectorField& QLes,
        volVectorField& QRans,
        volScalarField& mrEpsilonAvg,
        volSymmTensorField& mrRAvg,
        IOdictionary& relaxParameters
    );
      
      
      // Destructor
      ~relaxForcing();
      
      // Member Functions (public)
    
      //- Update all quantities and forcing terms
      void correct();

    // Access
      const Time &runTime() const { return runTime_;}
      const volVectorField& UL() { return UL_;}
      const volVectorField& UAvg() { return UAvg_;}
      const volSymmTensorField& RAvg() { return RAvg_;}
      const volScalarField& lesFlag() {return lesFlagL_;}
      bool couplingStatus() {return bool(couplingStatus_);}
      double rampQ() {return double(rampQ_);}  // for use in impact
      volVectorField& QSgs() {return QSgs_;} // for use in impact
      
      // Check
      static void printVersion();


      // Edit

      // Write

      // Member Operators
      
      // Friend Functions
      
      // Friend Operators

      // IOstream Operators

  };


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "relaxForcingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
