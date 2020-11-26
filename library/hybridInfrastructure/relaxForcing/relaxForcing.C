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
#include "wordRe.H"

#define PRINT_BOUNDS

namespace Foam {

  // * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

    //- Update consistent quantities
    void relaxForcing::updateConsistent()
    {
      // Averaging scheme for variable x:
      // Avg(x)[n+1] = existWt * Avg(x)[n] + instantWt * x[n+1]
      // existWt = 1/(1+dt/T) is the weight for the averaged field from previous step
      // instantWt = 1-existWt is the weight for the newly added instantaneous field
      // Ref: Meneveau et al. JFM. vol 319, pp 364; 
      // <<A lagrangian dynamic SGS model of turbulence>>
      
        dimensionedScalar existWt = 1.0/(1.0 + (runTime_.deltaT()/averagingTime_) );
        dimensionedScalar instantWt  = 1 - existWt;
        volSymmTensorField D = dev(symm(fvc::grad(UL_)));
        
        //- Update consistent velocity [Avg:1]
        UAvg_ = existWt * UAvg_ + instantWt * UL_;
      
        //- Update consistent dissipation rate [Avg:2]
        // epsilon := 2 * sgsModelRef_.nu() * magSqr(D) - (sgsModelRef_.B() && D);
        if(epsilonSgsBy_ == "dudx")
        {
            epsilonAvg_ = existWt * epsilonAvg_
                +  2.0 * sgsModelRef_.nuEff() * (magSqr(D) * instantWt);
        }
        else if(epsilonSgsBy_ == "kdelta")  // For impactFoam only
        {
            epsilonAvg_ = existWt * epsilonAvg_
                + ( 2.0 * sgsModelRef_.nu() * magSqr(D) + sgsModelRef_.epsilon() ) * instantWt;
        }
        else
        {
            epsilonAvg_ = existWt * epsilonAvg_
                +  2.0 * sgsModelRef_.nuEff() * (magSqr(D) * instantWt);
        }
        
        if( rmDivTauSgs_ )
        {
            //- Sgs modeling in the form of forcing in momentum equation (e.g. ADMRT forcing)
            QSgsAvg_ = existWt * QSgsAvg_ + instantWt * QSgs_;
        }
        
        //- Total kinetic stress [Avg:3]
        RResolvedAvg_ = existWt * RResolvedAvg_ + 
            instantWt * symm( ( UL_ - UAvg_ ) * (UL_ - UAvg_) );
        
        //- SGS turbulent stress [Avg:4]
        RSgsAvg_ = existWt * RSgsAvg_ + instantWt * sgsModelRef_.B();
        
        //- Average turbulent stress (needed by RANS model and force computation)
        RAvg_ = RResolvedAvg_ + RSgsAvg_;

        //- SGS turbulent energy [Avg:5] (need by resolution models, access via lookup)
        kSgsAvg_ = existWt * kSgsAvg_ + instantWt * sgsModelRef_.k();

        kAvg_ = 0.5 * tr(RAvg_); // only for output
        kLES_ = 0.5 * ( ( UL_ - UAvg_ ) & ( UL_ - UAvg_ ) ) + sgsModelRef_.k();

    } // End of updateConsistent


  //- Compute forcing terms QLES/QRans and rampQ
  void relaxForcing::computeForcing()
    {                
        volVectorField uprimeScale =  UL_ - UAvg_ ;
        
        // Velocity relax term + TKE/RR relax term
        // Note: mag(tr(.)) is not necessary. Just in case some thing goes
        // wrong during interpolation
        
        if (RConsisLevel_ == "full")
        {
            volSymmTensorField Gij = scalar(2.0) *
                ( meshTalk_.RR() - RAvg_ )
                /
                ( mag( tr(meshTalk_.RR()) ) + mag( tr(RAvg_) ) + ksmall_ );

            QLes_ = 
                (meshTalk_.UR() - UAvg_) / URelaxTime_ 
                + ( Gij & uprimeScale / RRelaxTime_ );
        }
        else if (RConsisLevel_ == "diag")
        {
            volSymmTensorField Gij = scalar(2.0) *
                ( meshTalk_.RR() - RAvg_ )
                /
                ( mag( tr(meshTalk_.RR()) ) + mag( tr(RAvg_) ) + ksmall_ );

            Gij.replace(symmTensor::XY, 0.0);
            Gij.replace(symmTensor::XZ, 0.0);
            Gij.replace(symmTensor::YZ, 0.0);

            QLes_ = 
                (meshTalk_.UR() - UAvg_) / URelaxTime_ 
                + ( Gij & uprimeScale / RRelaxTime_ );
            
        }
        else if (RConsisLevel_ == "shear")
        {
            volSymmTensorField Gij 
                (
                    IOobject
                    (
                        "relaxForcing::Gij",
                        runTime_.timeName(),
                        meshL_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    meshL_,
                    dimensionedSymmTensor("zero", dimless, symmTensor::zero)
                );
            
            Gij.replace
                (
                    symmTensor::XY,
                    scalar(2.0) *
                    (2.0* meshTalk_.RR().component(symmTensor::XY) - RAvg_.component(symmTensor::XY) )
                    /
                    ( mag( tr(meshTalk_.RR()) ) + mag( tr(RAvg_) ) + ksmall_ )
                );

            QLes_ = 
                (meshTalk_.UR() - UAvg_) / URelaxTime_ 
                + ( Gij & uprimeScale / RRelaxTime_ );
        }
        else if (RConsisLevel_ == "k")
        {
            volScalarField Gk = scalar(2.0) *
                ( 2.0 * meshTalk_.kR() - mag( tr(RAvg_) ) )
                /
                ( 2.0 * meshTalk_.kR() + mag( tr(RAvg_) ) + ksmall_ );
            
            QLes_ = 
                (meshTalk_.UR() - UAvg_) / URelaxTime_ 
                + ( Gk * uprimeScale / RRelaxTime_ );
        }
        else if (RConsisLevel_ == "fullDirect")
        {
            checkRCorrPtr();
            volSymmTensorField & Gij = RCorrPtr_();
            Gij =  meshTalk_.RR() - RAvg_;

            QLes_ = 
                // (meshTalk_.UR() - UAvg_) / URelaxTime_ 
                - fvc::div(Gij, "div(RR)");

            Gij = Gij * (scalar(1) - lesFlagL_);
        }
        else if (RConsisLevel_ == "shearDirect")
        {
            checkRCorrPtr();
            volSymmTensorField & Gij = RCorrPtr_();
            Gij =  meshTalk_.RR() - RAvg_;

            Gij.replace(symmTensor::XX, 0.0);
            Gij.replace(symmTensor::YY, 0.0);
            Gij.replace(symmTensor::ZZ, 0.0);
            Gij.replace(symmTensor::YZ, 0.0);
            Gij.replace(symmTensor::XZ, 0.0);

            QLes_ = 
                // (meshTalk_.UR() - UAvg_) / URelaxTime_ 
                - fvc::div(Gij, "div(RR)");

            Gij = Gij * (scalar(1) - lesFlagL_);
        }
        else if (RConsisLevel_ == "flexDirect")
        {
            checkRCorrPtr();
            volSymmTensorField & Gij = RCorrPtr_();
            Gij =  meshTalk_.RR() - RAvg_;

            if(flexMask_.xx() < SMALL)
            { 
                Gij.replace(symmTensor::XX, 0.0); 
            }
            
            if(flexMask_.yy() < SMALL)
            { 
                Gij.replace(symmTensor::YY, 0.0);
            }
            
            if(flexMask_.zz() < SMALL)
            { 
                Gij.replace(symmTensor::ZZ, 0.0);
            }
            
            if(flexMask_.xy() < SMALL)
            { 
                Gij.replace(symmTensor::XY, 0.0);
            }
            
            if(flexMask_.xz() < SMALL)
            { 
                Gij.replace(symmTensor::XZ, 0.0);
            }
            
            if(flexMask_.yz() < SMALL)
            { 
                Gij.replace(symmTensor::YZ, 0.0);
            }

            QLes_ = 
                // (meshTalk_.UR() - UAvg_) / URelaxTime_ 
                - fvc::div(Gij, "div(RR)");

            Gij = Gij * (scalar(1) - lesFlagL_);
        }
        else if (RConsisLevel_ == "diagDirect")
        {
            checkRCorrPtr();
            volSymmTensorField & Gij = RCorrPtr_();
            Gij =  meshTalk_.RR() - RAvg_;
            
            Gij.replace(symmTensor::XY, 0.0);
            Gij.replace(symmTensor::XZ, 0.0);
            Gij.replace(symmTensor::YZ, 0.0);
            
            QLes_ = 
                - fvc::div(Gij, "div(RR)");

            Gij = Gij * (scalar(1) - lesFlagL_);
        }
        else if (RConsisLevel_ == "laplacian") 
        {
            surfaceVectorField dudn ("dudn", fvc::snGrad(UAvg_ - meshTalk_.UR(), "DU"));
            surfaceScalarField rNeiDist("deltaCoeffs",  meshL_.surfaceInterpolation::deltaCoeffs());
            surfaceScalarField dFlag
                (
                    "dFlag", 
                    1.0 - pos(
                        mag(fvc::snGrad(lesFlagL_, "lesFlag") / rNeiDist)
                        - dimensionedScalar("modSmall", dimless, 0.1)
                    )
                );

            surfaceVectorField deltaFlux("deltaFlux",   nuRelaxFactorL_ * nu0_ * dudn * dFlag * meshL_.magSf());
            
            QLes_ = fvc::div(deltaFlux);
        }
        else
        {
            FatalErrorIn("relaxForcing::checkRestart()")
                << "RConsisLevel " << RConsisLevel_ << " not recognized." << nl
                << "Should be one of the following: " << endl
                << "k (default), full, shear, diag, laplacian" << endl
                << abort(FatalError);
        }
        
        if (RConsisLevel_ == "laplacian") 
        {
            surfaceVectorField dudn ("dudn", fvc::snGrad( UR_ - meshTalk_.UAvg(), "DU"));

            surfaceScalarField rNeiDist("deltaCoeffs",  meshR_.surfaceInterpolation::deltaCoeffs());
            surfaceScalarField dFlag
                (
                    "dFlag", 
                    1.0 - pos(
                        mag(fvc::snGrad(lesFlagR_, "lesFlag") / rNeiDist)
                        - dimensionedScalar("modSmall", dimless, 0.1)
                    )
                );

            surfaceVectorField deltaFlux("deltaFlux",  nuRelaxFactorR_ * nu0_ * dudn * dFlag * meshR_.magSf());
            QRans_ = fvc::div(deltaFlux);
        }
        else
        {
            QRans_ = ( meshTalk_.UAvg() - UR_ ) / MRelaxTime_;               
        }
        
        // Remove the Sgs forcing already in place, e.g. for ADMRT
        // And then mask the relax forcing fields
        if( rmDivTauSgs_ )
        {
            QLes_ = (QLes_ - QSgsAvg_) * (scalar(1) - lesFlagL_);
        }
        else
        {
            QLes_ = QLes_  * (scalar(1) - lesFlagL_);
        }

        QRans_ = QRans_ * lesFlagR_;
        
        if(neutralQMomentum_)
        {
            QLes_  -= QLes_.weightedAverage( meshL_.V() );
            QRans_ -= QRans_.weightedAverage( meshR_.V() );
        }
        
        {
            dimensionedVector QLAvg = QLes_.weightedAverage( meshL_.V() );
            dimensionedVector QRAvg = QRans_.weightedAverage( meshR_.V() ); 
            
            Info << "Pre diagnosis: " << nl 
                 << "    Avg QLes  (all) = " << QLAvg.value() << nl
                 << "    Avg QRans (all) = " << QRAvg.value() << endl;
        }
        
        // If during the transition period:
        // apply a linear ramp on Q
        if(runTime_.time() < (couplingStartTime_ + gradualCouplingDuration_))
        {
          rampQ_ = (
                    (
                     runTime_.time() - couplingStartTime_
                     ) / gradualCouplingDuration_
                    ).value();
          
          rampQ_ = Foam::min(Foam::max(rampQ_, scalar(0)), scalar(1));
          QLes_ = QLes_ * rampQ_;
          QRans_ = QRans_ * rampQ_;
        } 
      else
      { 
          rampQ_ = scalar(1);
      }
      
      consistentDict_.set("rampQ", rampQ_);
    }


  //- Auxillary functions
  
  //- Initialize consistent dictionary (coupling status and Ar)
  void  relaxForcing::initConsistentDict()
  {
    if (!consistentDict_.headerOk())
      {
        Info << "consistentDict not found." << endl;
        consistentDict_.add("couplingStatus", couplingStatus_);
        consistentDict_.add("rampQ", rampQ_);
      }
    else
      {
        Info << "\nReading couplingStatus from current directory. \n" << endl;
        consistentDict_.readIfPresent("couplingStatus", couplingStatus_); 
        consistentDict_.readIfPresent("rampQ", rampQ_);
      }       
    Info << "consistentDict initialized as follows (either from file or default):" << endl;
    Info << "consistentDict" << consistentDict_ << endl;
  }


  void relaxForcing::updateCouplingStatus()
  {
    if (
        (! couplingStatus_) &&
        LesRansCoupling_    &&
        runTime_.time() > (couplingStartTime_ - runTime_.deltaT() * 0.01)
        )
      {
        Info << "\n<< LES/RANS coupling turned on from Time = " 
             << runTime_.time().value() << " >>\n" << endl;
        couplingStatus_ = true;
        // Change dict entry.
        consistentDict_.set("couplingStatus", couplingStatus_);
      }
  }     
    

  //- Diagnosis of difference of UL & UR
  void relaxForcing::diagnosisDiffU()
  {

      if(consistencyDiagnosis_)
      {
          volScalarField diffULField = Foam::mag( UAvg_ - meshTalk_.UR() );
          
          dimensionedScalar diffUL = diffULField.weightedAverage( meshL_.V() * lesFlagL_ + vsmall_ ); 
          
          volScalarField diffURField = Foam::mag( meshTalk_.UAvg() - UR_ );
          
          volScalarField rasFieldR = 1.0 - lesFlagR_;
          dimensionedScalar diffUR = diffURField.weightedAverage( meshR_.V() * rasFieldR + vsmall_ ); 
          
          Info << "Consistency diagnosis: " << nl
               << "    RANS-LES velocity deviation : "
               << "LES zone diff(U) = " << diffUL.value() << ", "
               << "RAS zone diff(U) = " << diffUR.value() << endl;
          
          volScalarField resolvedTKE = scalar(0.5) * Foam::magSqr(UL_ - UAvg_);
          dimensionedScalar avgTKE = resolvedTKE.weightedAverage(meshL_.V() + vsmall_);
          
          Info << "    Resolved LES TKE: " << avgTKE.value() << endl;
      }

      if(forcingDiagnosis_)
      {
          // relaxation forcing
          volScalarField rasFieldL = 1.0 - lesFlagL_;
          dimensionedVector QLAvg = QLes_.weightedAverage( meshL_.V() * rasFieldL + vsmall_ );
          dimensionedVector QRAvg = QRans_.weightedAverage( meshR_.V() * lesFlagR_ + vsmall_ ); 
          
          Info << "Forcing diagnosis: " << nl 
               << "    Avg QLes in RANS zone = " << QLAvg.value() << nl
               << "    Avg QRans in LES zone = " << QRAvg.value() << endl;

          // shear stresses, LES
          const fvPatchList& patchesL = meshL_.boundary();
        
          forAll(patchesL, patchi)
          {
              const fvPatch& currPatch = patchesL[patchi];
              
              if (typeid(currPatch) == typeid(wallFvPatch))
              {
                  scalar area = gSum(meshL_.magSf().boundaryField()[patchi]);
                  vector sumField(vector::zero);
                  vector sumFieldNu(vector::zero);
                  
                  if (area > 0)
                  {
                      sumField = gSum
                          (
                              meshL_.magSf().boundaryField()[patchi]
                              * (sgsModelRef_.nuEff())().boundaryField()[patchi] 
                              * UL_.boundaryField()[patchi].snGrad()
                          ) / area;
                      
                      sumFieldNu = gSum
                          (
                              meshL_.magSf().boundaryField()[patchi]
                              * sgsModelRef_.nu().boundaryField()[patchi] 
                              * UL_.boundaryField()[patchi].snGrad()
                          ) / area;
                  }
                  Info << "    Avg Tauw(L) over patch "
                       <<  currPatch.name() << " = " << sumField 
                       <<  " viscous only = " << sumFieldNu << endl;
            }
          }

          // shear stresses, RANS
          const fvPatchList& patchesR = meshR_.boundary();
        
          forAll(patchesR, patchi)
          {
              const fvPatch& currPatch = patchesR[patchi];
              
              if (typeid(currPatch) == typeid(wallFvPatch))
              {
                  scalar area = gSum(meshR_.magSf().boundaryField()[patchi]);
                  vector sumField(vector::zero);
                  vector sumFieldNu(vector::zero);
                  
                  if (area > 0)
                  {
                      sumField = gSum
                          (
                              meshR_.magSf().boundaryField()[patchi]
                              * (ransModelRef_.nuEff())().boundaryField()[patchi] 
                              * UR_.boundaryField()[patchi].snGrad()
                          ) / area;

                      sumFieldNu = gSum
                          (
                              meshR_.magSf().boundaryField()[patchi]
                              * ransModelRef_.nu().boundaryField()[patchi] 
                              * UR_.boundaryField()[patchi].snGrad()
                          ) / area;
                  }
                  Info<< "    Avg Tauw(R) over patch "
                      <<  currPatch.name() << " = " << sumField       
                      <<  " viscous only = " << sumFieldNu << endl;              
              }
          }
          
      }
  }


void relaxForcing::checkRestart()
{
    Switch isRestart = false;
    if(runTime_.time() > runTime_.deltaT() * 0.5)
    { isRestart = true; }
    
    couplingDict_.readIfPresent("isRestart", isRestart);
    
    if(isRestart)
    {
        Info << "*** CONTINUATION RUN!" << nl << endl;
        
        int headerCount = 
            UAvg_.headerOk()
            + RResolvedAvg_.headerOk() + RSgsAvg_.headerOk() 
            + RAvg_.headerOk() + epsilonAvg_.headerOk() + kSgsAvg_.headerOk();
        
        if (headerCount < 6)
        {
            FatalErrorIn("relaxForcing::checkRestart()")
                << "Restart failed because only " << headerCount 
                << " field are found!" << nl
                << "Debug Info: "  << nl
                << "UAvg: " << UAvg_.headerOk()
                << "  RResolvedAvg: " << RResolvedAvg_.headerOk() 
                << "  RSgsAvg: " << RSgsAvg_.headerOk() 
                << "  RAvg: " << RAvg_.headerOk()
                << "  epsilonAvg: " << epsilonAvg_.headerOk()
                << "  kSgsAvg: " << kSgsAvg_.headerOk()
                << abort(FatalError);
        }
    }
    else // Re-initialize fields
    {
        Info << "*** FRESH SIMULATION!" << nl << endl;
        
        UAvg_ = UL_*1.0;
        
        RResolvedAvg_ *= 0.0;
        RSgsAvg_ *=  0.0;
        RAvg_ *= 0.0;
        
        epsilonAvg_ *= 0;
        kSgsAvg_ = sgsModelRef_.k() * 1.0;
    }
}


//- Print relaxation parameters
void relaxForcing::printRelaxParameters()
{
    // Bound the coupling time parameters
    if ( couplingStartTime_ < runTime_.deltaT() 
    || gradualCouplingDuration_ < runTime_.deltaT() )
    {
        FatalErrorIn("relaxForcing constructor: ") 
            << "Improperly specified time for coupling: "
                << "couplingStartTime = " << couplingStartTime_.value() << nl
                << "gradualCouplingDuation = " << gradualCouplingDuration_.value()
                << endl << abort(FatalError);
    }
    else
    {
        Info << "Summary of relaxation parameters: \n" 
             << "**************************************\n"
             << "timeScales" << timeScaleDict_ << "\n"
             << "couplingOptions" << couplingDict_
             << "**************************************\n"
             << endl;
    }
}


void relaxForcing::checkRCorrPtr()
{ 
    if(! RCorrPtr_.valid())
    {
        FatalErrorIn("relaxForcing::computeForcing()")
            << "RConsisLevel = " << RConsisLevel_  << nl
            << "but RCorrPtr_ is not valie: " << endl
            << abort(FatalError);
    }
}
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //  
#include "relaxForcingConstructor.H"

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

  relaxForcing::~relaxForcing()
  {}



  // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
  void relaxForcing::correct()
  {
      label timeIndex = runTime_.timeIndex();

      updateCouplingStatus(); // turn on coupling if necessary
      updateConsistent();    
      
      // update the consistent quantities on RANS grid
      if( timeIndex % min(mapR2LEvery_, mapL2REvery_) == 0 )
      {
          Info << "lesFlag updated." << endl;
          turbFlag_ -> evaluateFlag(); 
      }

      // Update R2L fields (needed for update QL)
      if( timeIndex % mapR2LEvery_ == 0)
      {
          Info << "Mapping fields (R2L)." << endl;
          meshTalk_.interpolate("R2L");
      }
      else
      {
          Info << "R2L mapping skipped: #"
               << timeIndex % mapR2LEvery_ << endl;
      }

      // Update L2R fields (needed for update QR)
      if( timeIndex % mapL2REvery_ == 0)
      {
          Info << "Mapping fields (L2R)." << endl;
          meshTalk_.interpolate("L2R");
      }
      else
      {
          Info << "L2R mapping skipped: #"
               << timeIndex % mapL2REvery_ << endl;
      }

      // update forcing (every time step)
      if ( LesRansCoupling_ && couplingStatus_ ) // coupling enabled and time ready
      {
          computeForcing(); // forcing on momentum equations
          ransModelRef_.enforceFields(turbRelaxTime_, rampQ_, enforceMode_);
      }
      else
      {
          QLes_ *= 0.0; QRans_ *= 0.0;
          Info << "LES/RANS coupling disabled." << endl;
      }

   
      // Display diagnosis information of UL/UR fields (difference)
      ransModelRef_.diagnosis();
      diagnosisDiffU();
  }

  // * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

  // * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

  // * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace Foam

// ************************************************************************* //
