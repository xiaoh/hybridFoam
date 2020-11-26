// Last update:
// Time-stamp: <2010-03-31 22:29:29 xiaoh>
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

  \*----------------------------------------------------------------------------*/

#include "relaxForcing.H"
#include "wallDist.H"


// #define PRINT_BOUNDS
#define DEBUG_Q

namespace Foam {

  // * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

    //- Update consistent quantities
    void relaxForcing::updateConsistent()
    {
      //- Ar of previous step
      scalar Ar0 = Ar; 
      // Accumulate ratio
      Ar = ratio * Ar0 + 1.0;
      consistentDict_.set("accumulatedRatio", Ar);

      // Update consistent velocity
      UAvg_ = (ratio * Ar0 * UAvg_ + UL_)/Ar;

      // Update consistent dissipation rate
      volSymmTensorField D = dev(symm(fvc::grad(UL_)));
      //       epsilonAvg_ = (ratio * Ar0 * epsilonAvg_ + 
      //                      2 * sgsModelRef_.nu() * magSqr(D)
      //                      - (sgsModelRef_.B() && D))/Ar;
      
      epsilonAvg_ = (ratio * Ar0 * epsilonAvg_ + 
                     2 * sgsModelRef_.nuEff() * magSqr(D))/Ar;
      
      if(preSpecifiedZones_) 
        {
          Info << "Skipping consistent turburlent quantity update."<<endl;
          return;
        }

      // Update consistent quantities
      //- Total kinetic stress
      RResolvedAvg_ = (ratio * Ar0 * RResolvedAvg_ + symm(UL_ * UL_))/Ar;

      //- SGS turbulent stress
      RSgsAvg_ = (ratio * Ar0 * RSgsAvg_ + sgsModelRef_.B())/Ar;

      //- Consistent turbulent stress
      RAvg_ = RResolvedAvg_ + RSgsAvg_ - symm(UAvg_ * UAvg_);

      //- SGS turbulent energy
      kSgsAvg_ = (ratio * Ar0 * kSgsAvg_
                  + sgsModelRef_.k())/Ar;

      #ifdef PRINT_BOUNDS
      Info << "@@@ sgs TKE by kSgsAvg: " 
           << Foam::gMin(kSgsAvg_).value() << ", " 
           << Foam::gMax(kSgsAvg_).value() << endl;
      Info << "@@@ range(k from RAvg): " 
           << Foam::gMin(0.5*tr(RAvg_)).value() << ", "
           << Foam::gMax(0.5*tr(RAvg_)).value()  << endl;
      #endif

    } // End of updateConsistent
  
  
  //- Update integral length scale (or other criteria) and flag
    void relaxForcing::evaluateResolution()
    {

      dimensionedScalar kzero("kzero", dimVelocity*dimVelocity, 0.0);
      dimensionedScalar epsSmall
        (
         "epsSmall", dimVelocity*dimVelocity/dimTime, SMALL*10.0
         );
      
      if(Foam::min(Foam::mag(epsilonAvg_)).value() < SMALL*10.0)
        Info << "*** WARNING: epsilonAvg is too small. "
             << "I will bound it with 10*SMALL." << endl;
      
      if(Foam::min(epsilonAvg_).value() < 0.0)
        Info << "*** WARNING: epsilonAvg is negative. " 
             << "I will bound it with 10*SMALL" << endl;
      
      // Evaluate according to length criterion
          
      // Total TKE computed in a "postive-definite" way
      // <UL * UL> - UAvg * UAvg + kSgsAvg,
      // where <*> is the consistent/LES operator
      volScalarField kPD = scalar(0.5) * (tr(RResolvedAvg_) - magSqr(UAvg_)) 
                           + kSgsAvg_;
      
      #ifdef PRINT_BOUNDS
      Info << "@@@ range kPD: " 
           << Foam::gMin(kPD).value() << ", "
           << Foam::gMax(kPD).value() << endl;
      Info << "@@@ resolved TKE: "
           << Foam::gMin(scalar(0.5) * 
                        (tr(RResolvedAvg_) - magSqr(UAvg_))).value()
           << ", "
           << Foam::gMax(scalar(0.5) * 
                        (tr(RResolvedAvg_) - magSqr(UAvg_))).value()
           << endl;
      #endif

      if(gMin(kPD) < 0.0)
        Info << "*** WARNING: TKE is negative, min = " 
             << gMin(kPD) << endl;
      
      kPD = Foam::max(kPD, kzero); // bound k by zero

      turbLength_ = 
        Foam::sqrt(pow3(kPD)) / Foam::max(epsilonAvg_, epsSmall); 

      /* //<@@@ 
      Info << "@@@ range turbLength " 
           << Foam::min(turbLength_).value() << ", "
           << Foam::max(turbLength_).value() << endl;
      //@@@> */

      resIndicatorL_.internalField() = 
        Foam::pow((Foam::pow3(turbLength_.internalField()) 
                   / (mesh_.V())), 0.333);

      // Update lesFlagL field with latency strategy
      lesFlagL_.internalField() += 
        (
         ( 
          pos(Foam::pow3(turbLength_.internalField()
                         /TLAdjustCoeff_.value()
                         *resLLower_.value())
              -   (mesh_.V())
              )
          * (1.0-lesFlagL_.internalField())
          )
         -
         (
          pos( 
              (mesh_.V()) -
              Foam::pow3(turbLength_.internalField()
                         /TLAdjustCoeff_.value()
                         *resLUpper_.value())
             )
          * lesFlagL_.internalField()
          )
         );
      
      // Evaluate according to energy criterion (Pope 2004): 
      // Resolve more than 80% of TKE
      dimensionedScalar kSmall("kSmall", dimVelocity*dimVelocity, scalar(SMALL*10.0) );
      resIndicatorK_ = 1.0 - Foam::max(kSmall, kSgsAvg_)
                    / Foam::max(kPD, Foam::max(kSmall, kSgsAvg_));
      
      //NOTE: kPD = scalar(0.5) * (tr(RResolvedAvg_) - magSqr(UAvg_)) 
      //                                   + kSgsAvg_;
      // Essentially, this bounds (tr(RResolvedAvg_) - magSqr(UAvg_)) by zero.

      // Update lesFlagK_ using latency strategy.
      lesFlagK_ += (
                    pos(resIndicatorK_ - resKUpper_.value()) * (1.0 - lesFlagK_) 
                    - pos(resKLower_.value() - resIndicatorK_) *      lesFlagK_ 
                    ); // Validated.
      
      // Calculate correlation
      volScalarField diffKL = Foam::mag(lesFlagK_ - lesFlagL_);
      scalar corrKL = scalar(1.0) - gAverage(diffKL);
      
      //        gAverage(Foam::sum(diffKL).value()/scalar(mesh_.nCells()));
      Info << "Correlation of L/K criteria: " << corrKL << endl;

      assertFlagInRange(lesFlagK_);
      assertFlagInRange(lesFlagL_);
    } // End of evaluateResolution
  

  //- Compute forcing terms QLES/QRans
  void relaxForcing::computeForcing()
    {
      if(LesRansCoupling_ && runTime_.time() > couplingStartTime_)
        { // Case 1: coupling on and has started.
          if (! couplingStatus_)
            { 
              Info << "\n<< LES/RANS coupling turned on from Time = " 
                   << runTime_.time().value() << " >>\n" << endl;
              couplingStatus_ = true;
              // Change dict entry.
              consistentDict_.set("couplingStatus", couplingStatus_);
            }
          
          QLes_ = (UR_ - UAvg_)/URelaxTime_ 
            + ( 
               ((RAvg_ - ransModelRef_.R())/ransModelRef_.k())
               & ((UAvg_ - UR_)/RRelaxTime_) 
               );
          QRans_ = (UAvg_ - UR_)/MRelaxTime_;
          
          // Mask relax forcing field
          QLes_ = QLes_ * (scalar(1) - *lesFlagPtr_);
          QRans_ = QRans_ * (*lesFlagPtr_);

          // Case 1a: During the transition period.
          // Apply a linear ramp on Q
          if(
             runTime_.time() <
             (couplingStartTime_ + gradualCouplingDuration_)
             )
            {
              rampQ_ = (
                       (
                        runTime_.time() - couplingStartTime_
                        ) / gradualCouplingDuration_
                       ).value();
              
              rampQ_ = Foam::min(Foam::max(rampQ_, scalar(0)), scalar(1));
              QLes_ = QLes_ * rampQ_;
              QRans_ = QRans_ * rampQ_;
              #ifdef DEBUG_Q
              Info << "@@@ time and rampQ: " << runTime_.time().value() 
                   << "  "  << rampQ_ << endl; 
              #endif
            } 
          else
            { rampQ_ = scalar(1);}
        }
      else
        { // Case 2: Coupling off, or not started yet.
          rampQ_ = scalar(0);
          QLes_ *= scalar(0);
          QRans_ *= scalar(0);
        }
    }


  //- Auxillary functions

  void relaxForcing::initConsistentFields()
  {
    int headerCount = UAvg_.headerOk()
      + RSgsAvg_.headerOk() 
      + RResolvedAvg_.headerOk() 
      + kSgsAvg_.headerOk()
      + epsilonAvg_.headerOk();
    
    Switch forcedInit(false); 

    if (headerCount == 0)
      {
        Info << "No consistent fields are found." << endl;
        Info << "Restarting consistent fields averaging." << endl;
        forcedInit = true;
      }
    else if(headerCount == 5)
      { 
        Info << "All 5 consistent fields loaded correctly!" << endl;
        Info << "Averaging continues from data saved previously at t = " 
             << runTime_.time().value() << " s \n" << endl;
        Info << "Computing RAvg from save data..." << endl;
        RAvg_ = RResolvedAvg_ + RSgsAvg_ - symm(UAvg_ * UAvg_);
        forcedInit = false;
      }
    else if( UAvg_.headerOk() 
             && epsilonAvg_.headerOk() 
             && preSpecifiedZones_ )
      {
        Info << "Pre-Specified LES/RANS zones, with UAvg & epsilonAvg loaded correctly!" << endl;
        Info << "Averaging continues from data saved previously at t = " 
             << runTime_.time().value() << " s \n" << endl;
        forcedInit = false;
      }
    else
      {
        Info << "Some consistent fields are missing. Healthy header count = "
             << headerCount << " (5 expected)." << endl;
        Info << "Cleaning all consistent fields."<< endl;
        UAvg_ = UL_;
        RSgsAvg_ *= 0.0;
        RResolvedAvg_ *= 0.0;
        kSgsAvg_ *= 0.0;
        epsilonAvg_ *= 0.0;
        lesFlagL_ = 1; lesFlagK_ = 1;
        forcedInit = true;
      }
    
    if (forcedInit)
      {
        Info << "Forced initialization of fields: ";
        Info << "RResolvedAvg, RSgsAvg, kSgsAvg and epsilonAvg.\n" << endl;
        // Initialize RResolvedAvg field
        RResolvedAvg_ = symm(UL_ * UL_);
        // The operation * 1.0 prevents aliasing.
        RSgsAvg_      = (sgsModelRef_.B() * 1.0);
        kSgsAvg_      = (sgsModelRef_.k() * 1.0);
        // Initialize epsilonAvg
        volSymmTensorField D = dev(symm(fvc::grad(UL_)));
        epsilonAvg_   = 2.0 * sgsModelRef_.nuEff() * magSqr(D);
      } 

    // Turn off consistent fields writing if so specified in Dictionary
    Switch writeConsistent
      (
       relaxParameters_.lookupOrAddDefault<Switch>("writeConsistentFields", true)
       );

    if(!writeConsistent)
      {
        Info << "*** Consistent fields will not be saved." << endl;
        RSgsAvg_.writeOpt() = IOobject::NO_WRITE; 
        RResolvedAvg_.writeOpt() = IOobject::NO_WRITE; 
        kSgsAvg_.writeOpt() = IOobject::NO_WRITE; 
        epsilonAvg_.writeOpt() = IOobject::NO_WRITE; 
      }
  }
  
  //- initialize consistent dictionary (coupling status and Ar)
  void  relaxForcing::initConsistentDict()
  {
    if (!consistentDict_.headerOk())
      {
        Info << "consistentDict not found." << endl;
        Info << "Initizlied couplingStatus with <false>." << endl
             << "Initizlied Ar with 1.\n" << endl;
        consistentDict_.add("couplingStatus", couplingStatus_);
        consistentDict_.add("accumulatedRatio", Ar);
      }
    else
      {
        Info << "Reading couplingStatus from current directory. \n" << endl;
        consistentDict_.readIfPresent("couplingStatus", couplingStatus_); 
        Info << "Coupling Status initizlied (either from file or defaults) as: <" 
             << couplingStatus_ << ">\n" << endl;
        consistentDict_.readIfPresent("accumulatedRatio", Ar); 
        Info << "Accumulated ratio initizlied (either from file or defaults) as: <" 
             << Ar << ">\n" << endl;
      }       
  }

  //- Assert UpperLimit > LowerLimit
  void relaxForcing::assertOrder(scalar& l, scalar& u)
  {
    if (l > u)
      {
        scalar tmp = u;
        u = l;
        l = tmp;
      }
  }

  //- Assert Flags are either 1 or 0
  void relaxForcing::assertFlagInRange(volScalarField& lesFlagX)
                        
  {
    scalar fmax = Foam::gMax(lesFlagX);
    scalar fmin = Foam::gMin(lesFlagX);
    word name =  lesFlagX.name();
    scalar small(1e-4);

    if ( 
        !(
          (mag(fmax-1.0) < small || mag(fmax) < small) &&
          (mag(fmin-1.0) < small || mag(fmin) < small) 
          )
       )
        FatalErrorIn("relaxForcing::assertFlagInRange()")
        << name << " should be either 0 or 1!" << endl
        << "Debug Info: " 
        << " max(" << name << ") = " << fmax << ", "
        << " min(" << name << ") = " << fmin << endl
        << abort(FatalError);
  }


  //- Choose resolution criterion
  void relaxForcing::chooseCriterion()
  {

    //- Potential inconsistency between here and createFieldsLR.H
    //- Need to be done in a more elegant way
    //- For now, the code need to be updated consistently anytime 
    //  wheevern a new resolution criteron is introduced.
    if(
       resolutionCriterion_ == "energy" || 
       resolutionCriterion_ == "Energy")
      {
        lesFlagPtr_ = & lesFlagK_;
      }
    else if 
      (
       resolutionCriterion_ == "length" || 
       resolutionCriterion_ == "Length"
       )
      {
        lesFlagPtr_ = & lesFlagL_;
      }
    else if
      (
       resolutionCriterion_ == "wallDist" || 
       resolutionCriterion_ == "walldist"
       )
      {
        lesFlagPtr_ = & lesFlagGeneral_;

        dimensionedScalar ransRangeTmp (couplingDict_.lookup("ransRange"));
         
        ransRange_ = ransRangeTmp;
        wallDist y(mesh_, true);
        lesFlagGeneral_.internalField() = pos(y.y() - ransRange_);
        
        preSpecifiedZones_ = true;
        switchFlagGeneralWriting();
      }
    else
      {
        lesFlagPtr_ = & lesFlagL_;
        Info << "Warning: relaxForcing::chooseCriterion:: unknown resolutionCriterion "
             << "in dictionary relaxationParameters.\n "
             << "Possible values: Length/Energy/wallDist"
             << "Default to Length criterion." << endl;
      }
    Info << "relaxForcing::chooseCriterion:: Resolution criterion used: "
         << lesFlagPtr_->name() << "\n" << endl;
  }


  //- Switch lesFlag writing options
  void relaxForcing:: switchFlagGeneralWriting()
  {
    Info << "Switching on lesFlagGeneral and off all coupling related fields."
         << endl;
    lesFlagK_.writeOpt() = IOobject::NO_WRITE; 
    lesFlagL_.writeOpt() = IOobject::NO_WRITE; 
    resIndicatorL_.writeOpt() = IOobject::NO_WRITE; 
    resIndicatorK_.writeOpt() = IOobject::NO_WRITE; 
    // Consistent Fields:
    RSgsAvg_.writeOpt() = IOobject::NO_WRITE; 
    RResolvedAvg_.writeOpt() = IOobject::NO_WRITE; 
    kSgsAvg_.writeOpt() = IOobject::NO_WRITE; 
    // epsilonAvg is needed for RANS computation.
    // epsilonAvg_.writeOpt() = IOobject::NO_WRITE; 
    
    // Switch on lesFlagGeneral
    lesFlagGeneral_.writeOpt() = IOobject::AUTO_WRITE; 
  }


  //- Diagnosis of difference of UL & UR
  void relaxForcing::diagnosisDiffU()
  {
    volScalarField diffULField =  
      Foam::mag((UAvg_ - UR_)* (*lesFlagPtr_));
    dimensionedScalar diffUL = gAverage(diffULField);
    
    volScalarField diffURField =
      Foam::mag((UAvg_ - UR_)*(scalar(1) - *lesFlagPtr_));
    dimensionedScalar diffUR = gAverage(diffURField); 

    Info << "RANS-LES velocity deviation : "
       << "LES zone diff(U) = " << diffUL.value() << ", "
       << "RAS zone diff(U) = " << diffUR.value() << endl;
  }


  //- Print relaxation parameters
  void relaxForcing::printRelaxParameters()
  {
    Info << "Summary of relaxation parameters: \n" 
         << "------------------------------------"
         << "------------------------------------\n"
         << "timeScales" << timeScaleDict_ << "\n"
         << "couplingOptions" << couplingDict_ << "\n"
         << "latencyCoeffs" << latencyDict_ << "\n"
         << "------------------------------------"
         << "------------------------------------\n"
         << endl;
  }
  
  // * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

  // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
  relaxForcing::relaxForcing(  
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
            IOdictionary& relaxParameters
			)
    :
    UL_(UL), UR_(UR),
    sgsModelRef_(sgsModelRef),
    ransModelRef_(ransModelRef),
    UAvg_(UAvg), RAvg_(RAvg), 
    RResolvedAvg_(RResolvedAvg), RSgsAvg_(RSgsAvg),
    epsilonAvg_(epsilonAvg),
    turbLength_(turbLength), 
    lesFlagL_(lesFlagL), lesFlagK_(lesFlagK), 
    lesFlagGeneral_(lesFlagGeneral),
    ransRange_("ransRange", dimLength, 0.0),
    resIndicatorL_(resIndicatorL), resIndicatorK_(resIndicatorK),
    QLes_(QLes), QRans_(QRans),
    runTime_(UL.time()),
    mesh_(UL.mesh()),
    // Dictionaries
    relaxParameters_(relaxParameters),
    timeScaleDict_(relaxParameters.subDictPtr("timeScales")),
    couplingDict_(relaxParameters.subDictPtr("couplingOptions")),
    latencyDict_(relaxParameters.subDictPtr("latencyCoeffs")),
    // Time scales:
    averagingTime_
        (
         dimensioned<scalar>::lookupOrAddToDict
         ("averagingTime", timeScaleDict_, 
          runTime_.endTime().value() * scalar(0.1),
          dimTime)
         ),

    URelaxTime_
        (
         dimensioned<scalar>::lookupOrAddToDict
         ("URelaxTime", timeScaleDict_, 
          averagingTime_.value() * scalar(0.2),
          dimTime)
         ),

    RRelaxTime_
        (
         dimensioned<scalar>::lookupOrAddToDict
         ("RRelaxTime", timeScaleDict_, 
          averagingTime_.value() * scalar(0.3),
          dimTime)
         ),

    MRelaxTime_
        (
         dimensioned<scalar>::lookupOrAddToDict
         ("MRelaxTime", timeScaleDict_, 
          averagingTime_.value() * scalar(0.2),
          dimTime)
         ),

    turbRelaxTime_
        (
         dimensioned<scalar>::lookupOrAddToDict
         ("turbRelaxTime", timeScaleDict_, 
          averagingTime_.value() * scalar(0.05),
          dimTime)
         ),

    // Coupling options:
    LesRansCoupling_
       (
        couplingDict_.lookupOrAddDefault<Switch>("LesRansCoupling", true)
        ),
    couplingStartTime_
        (
         dimensioned<scalar>::lookupOrAddToDict
         ("couplingStartTime", couplingDict_, 
          runTime_.endTime().value() * scalar(0.1),
          dimTime)
         ),

    gradualCouplingDuration_
        (
         dimensioned<scalar>::lookupOrAddToDict
         ("gradualCouplingDuration", couplingDict_,
          couplingStartTime_.value(), 
          dimTime)
         ),

    consistentDict_
            (
             IOobject
             (
              "consistentProperties",
              runTime_.timeName(),
              "consistent",
              mesh_,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE,
              true
              )
             ),

    rampQ_(1.0),
    enforceMode_
        (
         couplingDict_.lookupOrAddDefault<word>
         (
          "enforceMode", "both" 
          )
         ),

    resolutionCriterion_
        (
         couplingDict_.lookupOrAddDefault<word>
         (
          "resolutionCriterion", "length" 
          )
         ),
    TLAdjustCoeff_
        (
         dimensioned<scalar>::lookupOrAddToDict
         ("TLAdjustCoeff", couplingDict_, 12.0)
         ),
    resKLower_
        (
         dimensioned<scalar>::lookupOrAddToDict
         ("resKLower", latencyDict_, 0.78)
         ),
    resKUpper_
        (
         dimensioned<scalar>::lookupOrAddToDict
         ("resKUpper", latencyDict_, 0.82)
         ),
    resLLower_
        (
         dimensioned<scalar>::lookupOrAddToDict
         ("resLLower", latencyDict_, 0.9)
         ),
    resLUpper_
        (
         dimensioned<scalar>::lookupOrAddToDict
         ("resLUpper", latencyDict_, 1.1)
         ),

    couplingStatus_(false),
    preSpecifiedZones_(false),
    Ar(1.0),

    ratio(Foam::exp(-runTime_.deltaT()/averagingTime_).value()),
    kSgsAvg_
    (
     IOobject
     (
      "kSgsAvg",
      runTime_.timeName(),
      "consistent",
      mesh_,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
      ),
     sgsModelRef.k()
     )
  {
    // Point lesFlagPtr_ to the proper field;
    chooseCriterion();

    assertOrder(resKLower_.value(), resKUpper_.value());    
    assertOrder(resLLower_.value(), resLUpper_.value());

    initConsistentFields();

    // Bound the coupling time parameters
    couplingStartTime_ = Foam::max(couplingStartTime_, runTime_.deltaT());
    gradualCouplingDuration_ =
      Foam::max(gradualCouplingDuration_, runTime_.deltaT());

    initConsistentDict();
    printRelaxParameters();
  }

  
  
  // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

  relaxForcing::~relaxForcing()
  {}


  // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
  void relaxForcing::correct()
  {
    updateConsistent();

    if(preSpecifiedZones_)
      { Info << "Skipping resolution evaluation." << endl;}
    else
      { evaluateResolution();}

    computeForcing();
    
    if(LesRansCoupling_ && runTime_.time() > couplingStartTime_ ) 
      { ransModelRef_.enforceFields(turbRelaxTime_, rampQ_, enforceMode_); }

    // Display diagnosis information of UL/UR fields (difference)
    diagnosisDiffU();
  }


  // * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


  // * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


  // * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace Foam

// #include "relaxForcingIO.C"

// ************************************************************************* //
