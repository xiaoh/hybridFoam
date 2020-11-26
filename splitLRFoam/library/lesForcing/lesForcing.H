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
  lesForcing

  Description
  A class for the relaxation forcing in LES/RANS coupling

  SourceFiles
  lesForcing.H
  lesForcing.C


  \*---------------------------------------------------------------------------*/

#ifndef lesForcing_H
#define lesForcing_H

#include "Switch.H"
#include "muTurbulenceModel.H"
#include "IOdictionary.H"
#include "muLESModel.H"
#include "fvCFD.H"
#include "error.H"
#include "objectRegistry.H"
#include "Random.H"
#include "labelList.H"
#include "wallDist.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  /*---------------------------------------------------------------------------* \
    Class lesForcing Declaration
  \*---------------------------------------------------------------------------*/


  class lesForcing
  {
      
    // Private data
  private:

      //- References to velocity and pressure fields
      //  and turbulence models
      const volVectorField& U_;

      incompressible::muLESModel& sgsModelRef_;

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

      // References to the database and meshes
      const Time& runTime_;
      const fvMesh& mesh_;
      const label nCellsG_;
      
      // Dictionaries
      IOdictionary & relaxParameters_;
      dictionary     timeScaleDict_;
      dictionary     couplingDict_;

      // Time scale for averaging (consistent quantities)
      const dimensionedScalar averagingTime_;
      // Time scale for relaxation
      const dimensionedScalar URelaxTime_;
      const dimensionedScalar RRelaxTime_;

      // Whether to enforce LES fields on RAS, and
      // whether to add Q to the momentum equations
      const Switch LesRansCoupling_;
      word  RConsisLevel_;
      const dimensionedScalar couplingStartTime_;
      const dimensionedScalar gradualCouplingDuration_;
      IOdictionary consistentDict_;
      scalar rampQ_;
      const label mapR2LEvery_;
      Switch cleanFootprint_;

      Switch couplingStatus_;
      
      const dimensionedScalar ksmall_;
      const dimensionedScalar vsmall_;
      // Averaged SGS kinetic energy (come after mesh etc)
      volScalarField kSgsAvg_;

      volScalarField kAvg_;
      volScalarField kLES_; // total instantaneous LES energy, output only

      volScalarField lesFlag_;
      
      dimensionedScalar  ransRange_;
      volVectorField     mlUR_;
      volSymmTensorField mlRR_;
      volScalarField     mlkR_;

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

    //- Read from RANS fields (mapped on LES mesh).
    void updateRFields();
      
    //- Disallow default bitwise copy construct
      // lesForcing(const lesForcing&);

    //- Disallow default bitwise assignment
      // void operator=(const lesForcing&);


  public:


    // Constructors

    //- Construct from components
    lesForcing
    (            
        const volVectorField& U,
        Foam::incompressible::muLESModel& sgsModelRef,
        volVectorField& UAvg,
        volSymmTensorField& RAvg,
        volSymmTensorField& RResolvedAvg,
        volSymmTensorField& RSgsAvg,
        volScalarField& epsilonAvg,
        volVectorField& QLes,
        IOdictionary& relaxParameters
    );
      
      
      // Destructor
      ~lesForcing();
      
      // Member Functions (public)
      void barrier();
    
      //- Update all quantities and forcing terms
      void correct();

    // Access
      const Time &runTime() const { return runTime_;}
      const volVectorField& U() { return U_;}
      const volVectorField& UAvg() { return UAvg_;}
      const volSymmTensorField& RAvg() { return RAvg_;}
      
      // Check

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

// #include "lesForcingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
