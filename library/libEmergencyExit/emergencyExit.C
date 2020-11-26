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
  emergencyExit

  \*---------------------------------------------------------------------------*/

#include "emergencyExit.H"

namespace Foam {

  // * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

  
  // * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

  // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
  emergencyExit::emergencyExit(Time & runTime)
    :
    runTime_(runTime),
    
    emergencyExitDict_
    (
     IOobject
     (
      "emergencyExitDict",
      runTime_.constant(),
      runTime_,
      IOobject::READ_IF_PRESENT,
      IOobject::NO_WRITE
      )
     ),
  
    enabled_
       (
        emergencyExitDict_.lookupOrAddDefault<Switch>("emergencyExitSwitch", false)
        ),
    thresholdClockTime_
    (
     dimensioned<scalar>::lookupOrAddToDict
     ("thresholdClockTime", emergencyExitDict_, 
      27000, dimTime)
     ),
    forcedExitTime_
    (
     dimensioned<scalar>::lookupOrAddToDict
     ("forcedExitTime", emergencyExitDict_, 
      28000, dimTime)
     ),
    exitTimePattern_ 
    (emergencyExitDict_.lookupOrAddDefault<string>("exitTimePattern", "[1-9]+0")),
    wrPattern_(exitTimePattern_, wordRe::DETECT)
  {
    forcedExitTime_ = max(forcedExitTime_, thresholdClockTime_);
    Info << "emergencyExitDict: " <<  emergencyExitDict_ << endl;

    if (enabled_ && !Pstream::parRun())
    {
        deactivate();
        Info << "*** Emergency exit disabled as this is a serial run." << nl
            << endl;
    }
  }

  
  
  // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

  emergencyExit::~emergencyExit()
  {}


  // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool emergencyExit::execute(bool action)
  {
    bool alarmed = false;

    if (!enabled_) return alarmed;

    wordRe wrTimeName(runTime_.timeName());
    Switch patternMatched = wrPattern_.match(wrTimeName);
    
    if (
        runTime_.elapsedClockTime()  > thresholdClockTime_.value()
        && patternMatched
        )
      { 
        Info << "<<< thresholdClockTime reached and pattern matched. Envoke emergency exit at Time = " 
             << runTime_.timeName() << ", action = " << action << "  >>>" << endl;
        alarmed = true;

        if(action) { runTime_.writeAndEnd(); }

      }
    else if (runTime_.elapsedClockTime()  > forcedExitTime_.value())
      {
          Info << "<<< forcedExitTime reached. Envoke emergency exit at Time = " 
               << runTime_.timeName() << ", action = " << action << " >>>" << endl;
        alarmed = true;

        if(action) { runTime_.writeAndEnd(); }
      }
 
   return alarmed;
  }
  // * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


  // * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


  // * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace Foam


// ************************************************************************* //
