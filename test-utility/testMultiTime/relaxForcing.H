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
#include "LESModel.H"
#include "muTurbulenceModel.H"
#include "IOdictionary.H"
#include "muRASModel.H"
#include "fvCFD.H"
#include "error.H"
#include "objectRegistry.H"

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
    incompressible::LESModel& sgsModelRef_;
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
    // Integral length scale
    volScalarField& turbLength_;
    // Flag of LES/RANS domain
    volScalarField& lesFlagL_;
    volScalarField& lesFlagK_;
    volScalarField&  lesFlagGeneral_;
    dimensionedScalar ransRange_;
    //- Ref. to the actual flag used. Default to L
    volScalarField* lesFlagPtr_; 
    // Indicator of resolution (turbLength/meshSize)
    volScalarField& resIndicatorL_;
    volScalarField& resIndicatorK_;
    // Forcing for LES
    volVectorField& QLes_;
    // Forcing for RANS
    volVectorField& QRans_;

    // References to the database and meshes
    const Time& runTime_;
    const fvMesh& mesh_;
    
    // Dictionaries
    IOdictionary & relaxParameters_;
    dictionary           timeScaleDict_;
    dictionary           couplingDict_;
    dictionary           latencyDict_;

    // Time scale for averaging (consistent quantities)
    const dimensionedScalar averagingTime_;
    // Time scale for relaxation
    const dimensionedScalar URelaxTime_;
    const dimensionedScalar RRelaxTime_;
    const dimensionedScalar MRelaxTime_;
    const dimensionedScalar turbRelaxTime_;
    // Whether to enforce LES fields on RAS, and
    // whether to add Q to the momentum equations
    Switch LesRansCoupling_; 
    dimensionedScalar couplingStartTime_;
    dimensionedScalar gradualCouplingDuration_;
    IOdictionary consistentDict_;
    scalar rampQ_;
    word enforceMode_;

    // Use which criterion evaluate resolution
    word  resolutionCriterion_;
    const dimensionedScalar TLAdjustCoeff_;
    dimensionedScalar resKLower_;
    dimensionedScalar resKUpper_;
    dimensionedScalar resLLower_;
    dimensionedScalar resLUpper_;

    Switch couplingStatus_;
    Switch preSpecifiedZones_;

    // Accumulcated ratio 
    scalar Ar;
    // Ratio of geometric progression 
    // NOT run-time modifiable!
    const scalar ratio;

    // Averaged SGS kinetic energy
    volScalarField kSgsAvg_;
 
    // Private Member Functions

    //- Update integral length scale and flag
    void evaluateResolution();

    //- Update consistent quantities
    void updateConsistent();

    //- Compute forcing terms QLES/QRans
    void computeForcing();

    //- Check the consistent fields 
    //        and (re-)initialize them if necessary.
    void initConsistentFields();

    void initConsistentDict();

    //- Asserg lesFlagK and lesFlagL in the range [0 1]
    void assertFlagInRange(volScalarField & lesFlagX);
    
    //-assert that l < u. Switch if not.
    void assertOrder(scalar& l, scalar& u);
 
    // Switch field write options if wallDist criterion is chosen.
    void switchFlagGeneralWriting();

    //- Choose whether to use energy or length criterion
    void chooseCriterion();

    //- Diagnosis of difference of UL & UR
    void diagnosisDiffU();

    //- Print relaxation parameters
    void printRelaxParameters();

    //- Disallow default bitwise copy construct
    // relaxForcing(const relaxForcing&);

    //- Disallow default bitwise assignment
    // void operator=(const relaxForcing&);


  public:

    // Constructors

    //- Construct from components
    relaxForcing(            
			const volVectorField& UL,
			const volVectorField& UR,
            Foam::incompressible::LESModel& sgsModelRef,
            Foam::incompressible::muRASModel& ransModelRef,            
            volVectorField& UAvg,
            volSymmTensorField& RAvg,
            volSymmTensorField& RResolvedAvg,
            volSymmTensorField& RSgsAvg,
            volScalarField& epsilonAvg,
            volScalarField& turbLength,
            volScalarField& lesFlagL,
            volScalarField& lesFlagK,
            volScalarField& lesFlagGeneral,
            volScalarField& resIndicatorL,
            volScalarField& resIndicatorK,
            volVectorField& QLes,
            volVectorField& QRans,
            IOdictionary & relaxParameters
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

// #include "relaxForcingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
