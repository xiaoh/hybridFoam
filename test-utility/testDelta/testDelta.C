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

Application
    initSpecifiedField

Description
    Create fields in a specified way.

Note
   For purpose of testing utilities.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Random.H"
#include "pointFields.H"
#include "vectorList.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
  // #   include "readTransportProperties.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  typedef List<vectorList> vectorListList;

  pointField points = mesh.points();
  labelListList cellptLab = mesh.cellPoints();

  vectorListList* vertexPosPtr = new vectorListList(cellptLab.size());
  vectorListList & vertexPos = * vertexPosPtr;

  //  Info << "Points: " << points << endl;
  //  Info << "cellPoints: " << cellptLab << endl;

  forAll(cellptLab, cellI)
    {
      vertexPos[cellI].setSize(cellptLab[cellI].size()); 
      forAll(cellptLab[cellI], pI)
        {
          label p = cellptLab[cellI][pI];
          vertexPos[cellI][pI] = points[p];
        }
    }

  scalar tmpDist(-1);
  forAll(deltaL, cellI)
    {
      scalar maxDist = 0;
      forAll(vertexPos[cellI], pI) 
        {
          forAll(vertexPos[cellI], pJ)
            {
              if 
                ( pI!=pJ && 
                  (tmpDist=mag(vertexPos[cellI][pI]-vertexPos[cellI][pJ])) > maxDist 
                  )
                maxDist = tmpDist;
            }
        }
      
      deltaL.internalField()[cellI] = maxDist;
    }
    
    for(runTime++; !runTime.end(); runTime++)
    {
      //        Info<< "Time = " << runTime.timeName() << nl << endl;

        runTime.write();

    }


    return(0);
}


// ************************************************************************* //
