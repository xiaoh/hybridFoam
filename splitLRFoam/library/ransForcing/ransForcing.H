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
  ransForcing

  Description
  A class for the relaxation forcing in LES/RANS coupling

  SourceFiles
  ransForcing.H
  ransForcing.C


  \*---------------------------------------------------------------------------*/

#ifndef ransForcing_H
#define ransForcing_H

#include "Switch.H"
#include "muTurbulenceModel.H"
#include "IOdictionary.H"
#include "muRASModel.H"
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
    Class ransForcing Declaration
  \*---------------------------------------------------------------------------*/


  class ransForcing
  {
      
    // Private data
  private:

      //- References to velocity and pressure fields
      //  and turbulence models
      const volVectorField& U_;

      incompressible::muRASModel& ransModelRef_;

      // Forcing for RANS
      volVectorField& QRans_;

      // References to the database and meshes
      const Time& runTime_;
      const fvMesh& mesh_;
      const label nCellsG_;
      
      // Dictionaries
      IOdictionary & relaxParameters_;
      dictionary           timeScaleDict_;
      dictionary           couplingDict_;

      // Time scale for relaxation
      const dimensionedScalar MRelaxTime_;
      const dimensionedScalar turbRelaxTime_;
      // Whether to enforce LES fields on RAS, and
      // whether to add Q to the momentum equations
      const Switch LesRansCoupling_;

      const dimensionedScalar couplingStartTime_;
      const dimensionedScalar gradualCouplingDuration_;
      IOdictionary consistentDict_;
      scalar rampQ_;
      word enforceMode_;
      const label mapL2REvery_;
      Switch cleanFootprint_;
      
      Switch couplingStatus_;
      
      const dimensionedScalar ksmall_;
      const dimensionedScalar vsmall_;

      volScalarField lesFlag_;
      dimensionedScalar  ransRange_;
      volVectorField mrUAvg_;
      volScalarField & mrEpsilonAvg_;
      volSymmTensorField & mrRAvg_;

    // Private Member Functions

    //- Update integral length scale and flag
    void evaluateResolution();

    //- Compute forcing terms QLES/QRans
    void computeForcing();

    void initConsistentDict();

    //check fields are presents if restart; re-initialize fields if fresh start.
    void checkRestart();

    void updateCouplingStatus();

    void updateLFields();

    //- Diagnosis of difference of UL & UR
    void diagnosisDiffU();

    //- Print relaxation parameters
    void printRelaxParameters();

      
    //- Disallow default bitwise copy construct
      // ransForcing(const ransForcing&);

    //- Disallow default bitwise assignment
      // void operator=(const ransForcing&);


  public:


    // Constructors

    //- Construct from components
    ransForcing
    (            
        const volVectorField& UR,
        Foam::incompressible::muRASModel& ransModelRef,            
        volVectorField& QRans,
        volScalarField& mrEpsilonAvg,
        volSymmTensorField& mrRAvg,
        IOdictionary& relaxParameters
    );
      
      
      // Destructor
      ~ransForcing();
      
      // Member Functions (public)
      void barrier();    
      //- Update all quantities and forcing terms
      void correct();

    // Access
      const Time &runTime() const { return runTime_;}
      
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

// #include "ransForcingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
