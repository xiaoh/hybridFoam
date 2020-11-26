/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

Application
    postAverage

Eelco van Vliet

Description
    Post-processes data from flow calculations
    For each time: calculates the time average of a sequence of fields and 
    writes time time average in the directory



\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "fieldAverageFunctionObject.H"
//#include "dictionary.H"
//#include "IFstream.H"
//#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);
    runTime.setTime(timeDirs[0], 0);
#   include "createMesh.H"

   // what is required to initilise only the field in fieldAverage ?
// dictionary dict(IFstream("controlDict")());
// fieldAverageFunctionObject avefield() ;
// fieldAverage bla avefield.read() ;
// avefield.fieldAverage(dict);

    forAll(timeDirs, timeI)
    {
       runTime.setTime(timeDirs[timeI], timeI);
       Info<< "Adding fields for time " << runTime.timeName() << endl;
#      include "createFields.H"

      // Average fields over channel down to a line
      runTime.functionObjects().execute();
    }

    Info<< "\nEnd" << endl;


    return 0;
}

// ************************************************************************* //
