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

  Description
  A class for finishing current simulation and write fields when time is up.

  SourceFiles
  emergencyExit.H
  emergencyExit.C
  \*---------------------------------------------------------------------------*/


#ifndef   emergencyExit_H
#define   emergencyExit_H

#include "Switch.H"
#include "IOdictionary.H"
#include "fvCFD.H"
#include "string.H"
#include "wordRe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  /*---------------------------------------------------------------------------* \
    Class emergencyExit Declaration
  \*---------------------------------------------------------------------------*/

  class emergencyExit
  {
    // Private data
  private:
    Time & runTime_;
    IOdictionary emergencyExitDict_;
    Switch enabled_;
    dimensionedScalar thresholdClockTime_;
    dimensionedScalar forcedExitTime_;
    string  exitTimePattern_;
    wordRe wrPattern_;
 
    // Private Member Functions


  public:

    // Constructors

    //- Construct from components
    emergencyExit(Time& runTime);


    // Destructor
    ~emergencyExit();

    // Member Functions (public)
    
      //bool execute();
    bool execute( bool action=true );
    void activate() {enabled_ = true;}
    void deactivate() {enabled_ = false;}

    // Access

    
    // Check
    Switch enabled() { return enabled_;}
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

#endif

// ************************************************************************* //
