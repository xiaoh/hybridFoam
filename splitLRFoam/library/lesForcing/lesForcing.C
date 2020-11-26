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

  \*----------------------------------------------------------------------------*/



#include "lesForcing.H"
#include "IFstream.H"
#include "dictionary.H"
#include <unistd.h>

namespace Foam {

  // * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

    //- Update consistent quantities
    void lesForcing::updateConsistent()
    {
      // Averaging scheme for variable x:
      // Avg(x)[n+1] = existWt * Avg(x)[n] + instantWt * x[n+1]
      // existWt = 1/(1+dt/T) is the weight for the averaged field from previous step
      // instantWt = 1-existWt is the weight for the newly added instantaneous field
      // Ref: Meneveau et al. JFM. vol 319, pp 364; 
      // <<A lagrangian dynamic SGS model of turbulence>>
      
        dimensionedScalar existWt = 1.0/(1.0 + (runTime_.deltaT()/averagingTime_) );
        dimensionedScalar instantWt  = 1 - existWt;
        volSymmTensorField D = dev(symm(fvc::grad(U_)));
        
        //- Update consistent velocity [Avg:1]
        UAvg_ = existWt * UAvg_ + instantWt * U_;
      
        //- Update consistent dissipation rate [Avg:2]
        // epsilon := 2 * sgsModelRef_.nu() * magSqr(D) - (sgsModelRef_.B() && D);
        epsilonAvg_ = existWt * epsilonAvg_ 
            +  2.0 * sgsModelRef_.nuEff() * (magSqr(D) * instantWt);
        //  + ( 2.0 * sgsModelRef_.nu() * magSqr(D) + sgsModelRef_.epsilon() ) * instantWt;
        
        //- Total kinetic stress [Avg:3]
        RResolvedAvg_ = existWt * RResolvedAvg_ + 
            instantWt * symm( ( U_ - UAvg_ ) * (U_ - UAvg_) );
        
        //- SGS turbulent stress [Avg:4]
        RSgsAvg_ = existWt * RSgsAvg_ + instantWt * sgsModelRef_.B();
        
        //- Average turbulent stress (needed by RANS model and force computation)
        RAvg_ = RResolvedAvg_ + RSgsAvg_;

        //- SGS turbulent energy [Avg:5] (need by resolution models, access via lookup)
        kSgsAvg_ = existWt * kSgsAvg_ + instantWt * sgsModelRef_.k();

        kAvg_ = 0.5 * tr(RAvg_); // only for output
        kLES_ = 0.5 * ( ( U_ - UAvg_ ) & ( U_ - UAvg_ ) ) + sgsModelRef_.k();

    } // End of updateConsistent


  //- Compute forcing terms QLES/QRans and rampQ
  void lesForcing::computeForcing()
    {                
        volVectorField uprimeScale =  U_ - UAvg_ ;
        
        // Velocity relax term + TKE/RR relax term
        // Note: mag(tr(.)) is not necessary. Just in case some thing goes
        // wrong during interpolation
        
        if (RConsisLevel_ == "full")
        {
            volSymmTensorField Gij = scalar(2.0) *
                ( mlRR_ - RAvg_ )
                /
                ( mag( tr(mlRR_) ) + mag( tr(RAvg_) ) + ksmall_ );

            QLes_ = 
                (mlUR_ - UAvg_) / URelaxTime_ 
                + ( Gij & uprimeScale / RRelaxTime_ );
        }
        else if (RConsisLevel_ == "diag")
        {
            volSymmTensorField Gij = scalar(2.0) *
                ( mlRR_ - RAvg_ )
                /
                ( mag( tr(mlRR_) ) + mag( tr(RAvg_) ) + ksmall_ );

            Gij.replace(symmTensor::XY, 0.0);
            Gij.replace(symmTensor::XZ, 0.0);
            Gij.replace(symmTensor::YZ, 0.0);

            QLes_ = 
                (mlUR_ - UAvg_) / URelaxTime_ 
                + ( Gij & uprimeScale / RRelaxTime_ );
            
        }
        else if (RConsisLevel_ == "shear")
        {
            volSymmTensorField Gij 
                (
                    IOobject
                    (
                        "lesForcing::Gij",
                        runTime_.timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedSymmTensor("zero", dimless, symmTensor::zero)
                );
            
            Gij.replace
                (
                    symmTensor::XY,
                    scalar(2.0) *
                    (2.0* mlRR_.component(symmTensor::XY) - RAvg_.component(symmTensor::XY) )
                    /
                    ( mag( tr(mlRR_) ) + mag( tr(RAvg_) ) + ksmall_ )
                );

            QLes_ = 
                (mlUR_ - UAvg_) / URelaxTime_ 
                + ( Gij & uprimeScale / RRelaxTime_ );
        }
        else
        {
            volScalarField Gk = scalar(2.0) *
                ( 2.0 * mlkR_ - mag( tr(RAvg_) ) )
                /
                ( 2.0 * mlkR_ + mag( tr(RAvg_) ) + ksmall_ );
            
            QLes_ = 
                (mlUR_ - UAvg_) / URelaxTime_ 
                + ( Gk * uprimeScale / RRelaxTime_ );
        }
        
        // Mask relax forcing field
        QLes_ = QLes_ * (scalar(1) - lesFlag_);
      
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
        } 
      else
      { 
          rampQ_ = scalar(1);
      }
      
      consistentDict_.set("rampQ", rampQ_);
    }


  //- Auxillary functions
  
  //- Initialize consistent dictionary (coupling status and Ar)
  void  lesForcing::initConsistentDict()
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


  void lesForcing::updateCouplingStatus()
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
  void lesForcing::diagnosisDiffU()
  {
      volScalarField diffULField = Foam::mag( UAvg_ - mlUR_ );

      dimensionedScalar diffUL = diffULField.weightedAverage( mesh_.V() * lesFlag_ + vsmall_ ); 
      
      Info << "RANS-LES velocity deviation : "
           << "LES zone diff(U) = " << diffUL.value() << endl;

      
      volScalarField resolvedTKE = scalar(0.5) * Foam::magSqr(U_ - UAvg_);
      dimensionedScalar avgTKE = resolvedTKE.weightedAverage(mesh_.V() + vsmall_);
      
      Info << "Resolved LES TKE: " << avgTKE.value() << endl;
  }


void lesForcing::checkRestart()
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
            FatalErrorIn("lesForcing::checkRestart()")
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
        
        UAvg_ = U_*1.0;
        
        RResolvedAvg_ *= 0.0;
        RSgsAvg_ *=  0.0;
        RAvg_ *= 0.0;
        
        epsilonAvg_ *= 0;
        kSgsAvg_ = sgsModelRef_.k() * 1.0;
    }
}


//- Print relaxation parameters
void lesForcing::printRelaxParameters()
{
    // Bound the coupling time parameters
    if ( couplingStartTime_ < runTime_.deltaT() 
    || gradualCouplingDuration_ < runTime_.deltaT() )
    {
        FatalErrorIn("lesForcing constructor: ") 
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

void lesForcing::updateRFields() 
{
    mlUR_.instance() = runTime_.timeName();
    mlRR_.instance() = runTime_.timeName();
    if(mlUR_.headerOk())
    {
        Info << "Reading mlUR ..." << endl;
        IFstream str(mlUR_.filePath()); 
        dictionary fieldDict(str); 
        mlUR_.dimensionedInternalField().readField(fieldDict, "internalField");
    }
    else
    {
        FatalErrorIn("lesForcing::updateRFields()")
            << " mlUR is not found @ " << mlUR_.path()
            << abort(FatalError);
    }
    
    if(mlRR_.headerOk())
    {
        Info << "Reading mlRR ..." << endl;
        
        IFstream str(mlRR_.filePath()); 
        dictionary fieldDict(str); 
        mlRR_.dimensionedInternalField().readField(fieldDict, "internalField");
    }
    else
    {
        FatalErrorIn("lesForcing::updateRFields()")
            << " mlRR is not found @ " << mlRR_.path()
            << abort(FatalError);
    }

    mlkR_ = 0.5*tr(mlRR_);
    
    // Remove the time directory / fields written for communication
    if((!runTime_.outputTime()) && cleanFootprint_ && Pstream::parRun() ) 
    {
        Foam::rmDir(mlUR_.path());
    }
}
 
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //  
#include "lesForcingConstructor.H"

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

  lesForcing::~lesForcing()
  {}



  // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void lesForcing::barrier()
{
    // parallel barrier
    if(Pstream::parRun())
    { 
        Info << "=== Parallel barrier ===" << endl;
        
        if(Pstream::master())
        {
            label myno=Pstream::myProcNo();              
            for(label i = 1; i<Pstream::nProcs(); i++)
            {
                OPstream vStream(Pstream::blocking, i);
                vStream << myno;
            }
        }
        else
        {
            label mano=1000;
            IPstream vecStream(Pstream::blocking, 0);
            vecStream >> mano;
        }
    }
    else
    {
        return;
    }
}


  void lesForcing::correct()
  {
      label timeIndex = runTime_.timeIndex();

      updateCouplingStatus(); // turn on coupling if necessary
      updateConsistent();    

      fileName basePath;
      if(Pstream::parRun())
      {
          basePath = runTime_.path()/"../";
      }
      else
      {
          basePath = runTime_.path();
      }

      // Update R2L fields
      if( timeIndex % mapR2LEvery_ == 0)
      {
          Info << "Mapping fields (R2L)." << endl;
          
          UAvg_.write();
          RAvg_.write();
          epsilonAvg_.write();
          barrier();

          if(Pstream::master())
          {
              fileName hostReady = basePath / "hostReady" + runTime_.time().timeName();
              fileName commReady = basePath / "commReady" + runTime_.time().timeName();

              if(cleanFootprint_)
              {
                  scalar prevt = runTime_.time().value() - runTime_.time().deltaT().value() * mapR2LEvery_;
                  fileName hostPrev(basePath/"hostReady"+runTime_.time().timeName(prevt));
                  if(isFile(hostPrev)) { Foam::rm(hostPrev); }
              }
              
              system("touch " + hostReady);
              
              Info << "Waiting for communicator(L) ";
              while(1)
              {
                  if(isFile(commReady))
                  {
                      Info <<" Done! Update mlXR fields."  << endl;
                      break;
                  }
                  else
                  {
                      Info << ".";
                      Foam::sleep(2); // ::usleep(200000);
                  }
              }

          }

          barrier();
          updateRFields();
      }
      else
      {
          Info << "R2L mapping skipped: #"
               << timeIndex % mapR2LEvery_ << endl;
      }

      // update forcing (every time step)
      if ( LesRansCoupling_ && couplingStatus_ ) // coupling enabled and time ready
      {
          computeForcing(); // forcing on momentum equations
      }
      else
      {
          QLes_ *= 0.0;
          Info << "LES/RANS coupling disabled." << endl;
      }

   
      // Display diagnosis information of UL/UR fields (difference)
      diagnosisDiffU();
  }

  // * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

  // * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

  // * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace Foam

// ************************************************************************* //
