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

  \*----------------------------------------------------------------------------*/


#include "ransForcing.H"
#include "IFstream.H"
#include <unistd.h>

#define PRINT_BOUNDS
// #define DEBUG_Q

namespace Foam {

  // * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

  //- Compute forcing terms QRans and rampQ
  void ransForcing::computeForcing()
    {                
        QRans_ = ( mrUAvg_ - U_ ) / MRelaxTime_;
        
        // Mask relax forcing field
        QRans_ = QRans_ * lesFlag_;
      
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
  void  ransForcing::initConsistentDict()
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


  void ransForcing::updateCouplingStatus()
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
  void ransForcing::diagnosisDiffU()
  {
      volScalarField diffURField = Foam::mag( mrUAvg_ - U_ );

      volScalarField rasFieldR = 1.0 - lesFlag_;
      dimensionedScalar diffUR = diffURField.weightedAverage( mesh_.V() * rasFieldR + vsmall_ ); 
      
      Info << "RANS-LES velocity deviation : "
           << "RAS zone diff(U) = " << diffUR.value() << endl;
  }


void ransForcing::updateLFields() 
{
    mrUAvg_.instance() = runTime_.timeName();
    mrRAvg_.instance() = runTime_.timeName();
    mrEpsilonAvg_.instance() = runTime_.timeName();

    if(mrUAvg_.headerOk())
    {
        Info << "Reading mrUAvg ..." << endl;
        IFstream str(mrUAvg_.filePath()); 
        dictionary fieldDict(str); 
        mrUAvg_.dimensionedInternalField().readField(fieldDict, "internalField");
    }
    else
    {
        FatalErrorIn("ransForcing::updateLFields()")
            << " mrUAvg is not found @ " << mrUAvg_.path()
            << abort(FatalError);
    }

    if(mrRAvg_.headerOk())
    {
        Info << "Reading mrRAvg ..." << endl;
        IFstream str(mrRAvg_.filePath()); 
        dictionary fieldDict(str); 
        mrRAvg_.dimensionedInternalField().readField(fieldDict, "internalField");
    }
    else
    {
        FatalErrorIn("ransForcing::updateLFields()")
            << " mrRAvg is not found @ " << mrRAvg_.path()
            << abort(FatalError);
    }

    if(mrEpsilonAvg_.headerOk())
    {
        Info << "Reading mrEpsilonAvg ..." << endl;
        IFstream str(mrEpsilonAvg_.filePath()); 
        dictionary fieldDict(str); 
        mrEpsilonAvg_.dimensionedInternalField().readField(fieldDict, "internalField");
    }
    else
    {
        FatalErrorIn("ransForcing::updateLFields()")
            << " mrEpsilonAvg is not found @ " << mrEpsilonAvg_.path()
            << abort(FatalError);
    }

    if((!runTime_.outputTime()) && cleanFootprint_ && Pstream::parRun() )
    { 
        Foam::rmDir(mrUAvg_.path()); 
    }
}


void ransForcing::checkRestart()
{
    Switch isRestart = false;
    if(runTime_.time() > runTime_.deltaT() * 0.5)
    { isRestart = true; }
    
    couplingDict_.readIfPresent("isRestart", isRestart);
    
    if(isRestart)
    {
        Info << "*** CONTINUATION RUN!" << nl << endl;
    }
    else // Re-initialize fields
    {
        Info << "*** FRESH SIMULATION!" << nl << endl;
    }
}


//- Print relaxation parameters
void ransForcing::printRelaxParameters()
{
    // Bound the coupling time parameters
    if ( couplingStartTime_ < runTime_.deltaT() 
    || gradualCouplingDuration_ < runTime_.deltaT() )
    {
        FatalErrorIn("ransForcing constructor: ") 
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

 
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //  
#include "ransForcingConstructor.H"

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

  ransForcing::~ransForcing()
  {}



  // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void ransForcing::barrier()
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


  void ransForcing::correct()
  {
      label timeIndex = runTime_.timeIndex();

      updateCouplingStatus(); // turn on coupling if necessary

      fileName basePath;
      if(Pstream::parRun())
      {
          basePath = runTime_.path()/"../";

      }
      else
      {
          basePath = runTime_.path();
      }

      // Update L2R fields
      if( timeIndex % mapL2REvery_ == 0)
      {
          Info << "Mapping fields (L2R)." << endl;
          
          const volScalarField& epsilonR_(mesh_.lookupObject<const volScalarField>("epsilonR") );;
          const volSymmTensorField& RR_(mesh_.lookupObject<const volSymmTensorField>("RR") );
      
          U_.write();
          RR_.write();
          epsilonR_.write();
          barrier();

          if(Pstream::master())
          {
              fileName hostReady = basePath / "hostReady" + runTime_.time().timeName();
              fileName commReady = basePath / "commReady" + runTime_.time().timeName();

              if(cleanFootprint_)
              {
                  scalar prevt = runTime_.time().value() - runTime_.time().deltaT().value() * mapL2REvery_;
                  fileName hostPrev(basePath/"hostReady"+runTime_.time().timeName(prevt));
                  if(isFile(hostPrev)) { Foam::rm(hostPrev); }
              }

              system("touch " + hostReady);

              // wait for communicator to be ready
              Info << "Waiting for communicator(R) ";
              while(1)
              {
                  if(isFile(commReady))
                  {
                      Info <<" Done! Update mrXAvg fields."  << endl;
                      break;
                  }
                  else
                  {
                      Info << ".";
                      Foam::sleep(2);
                  }
              }
          }
          
          barrier();
          updateLFields();
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
          QRans_ *= 0.0;
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
