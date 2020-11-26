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
    ransCellToProcX

Description
    Make the cell to processor mapping file for RANS according to LES decomposition
    The info about LES decomposition is read from 0/cellDist field.


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "IFstream.H"
#include "OFstream.H"
#include "Switch.H"
#include "DynamicList.H"
#include "SortableList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("tolerateDuplicate", "");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMeshLR.H"
    #include "createFields.H"
    
    bool removeDuplicate  = ! args.optionFound("tolerateDuplicate");
    
    const  volVectorField & cellCentersR = meshR.C();
    
    List<boundBox> bbsMeshL(nProcs);
  
    for (int procI=0; procI<nProcs; procI++)
      {
        Info<< nl << "Reading LES mesh on processor " << procI << endl;

            Time runTimeLES
            (
                Time::controlDictName,
                runTime.path(),
                fileName(word("processor") + name(procI))
            );

            fvMesh meshSource
            (
                IOobject 
               (
                    fvMesh::defaultRegion,
                    runTimeLES.timeName(),
                    runTimeLES
                )
            );

            Info<< "mesh size: " << meshSource.nCells() << endl;

            bbsMeshL[procI] = meshSource.bounds();
            Info << "bound box for proc " << procI << ": " << meshSource.bounds() << endl;
        }

  // Cell ID on the overlap region of two bound boxes.
  DynamicList <label > borderCellIDs(meshR.nCells()/100+1);
  
  for (int procI=0; procI<nProcs; procI++)
    {
      forAll(cellDistR, celli)
        {
          const point & pt = cellCentersR[celli];
          if  (bbsMeshL[procI].contains(pt))
            {
              if (cellDistR[celli] < 0) // This cell has not been assigned
                cellDistR[celli] = procI; 
              else // already assigned to another processor
                borderCellIDs.append(celli);
                }
        }
    }


  if (removeDuplicate)
    {  // check cells assigned to multiple processors
      borderCellIDs.shrink();
  
      // Sort the border cell IDs and eliminate duplicate cell labels
      SortableList <label > sortedBorderCell(borderCellIDs);
      borderCellIDs.clear();
      label prevCell = -1;
      forAll(sortedBorderCell, cI)
        {
          label currCell = sortedBorderCell[cI];
          if( currCell != prevCell) borderCellIDs.append(currCell);
          prevCell = currCell;
        }
      Info << "\nNumber of (unique) cells on processor borders: " << borderCellIDs.size() << endl;
    }
  else
    Info << "\nNumber of cells (maybe with duplicate) on processor borders: " << borderCellIDs.size() << endl;


  // Assign the cells on the border according to strict criteria
  forAll(borderCellIDs, cI)
    {
      // Find the ID of cell on L mesh closet to this R cell
      label celliR = borderCellIDs[cI];
      procBorderFlag[celliR] = 1;
      label nearCi = meshL.findNearestCell(cellCentersR[celliR]);
      // Find the processor ID of this cell on L mesh
      cellDistR[celliR] = cellDist[nearCi];
    }
  
  // check unassigned cells
  if(Foam::min(cellDistR.internalField()) < 0)
    {
      cellDistR.write();
      Info << cellDistR.internalField() << endl;
      FatalErrorIn("ransCellToProc")
        << "Some cells in RANS mesh are not assigned to any processor." << abort(FatalError);
    }

  forAll(cellDistR, celli)
    {
      cellDecompositionR[celli] = cellDistR[celli];
    }
  
  Info<< nl << "Wrote decomposition files to " << nl
      << cellDecompositionR.objectPath() << nl
      << " for use in manual decomposition." << endl;
  
  cellDecompositionR.write();
  cellDistR.write();
  procBorderFlag.write();
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  
  Info<< nl << "cellDistR and procBorderFlag scalar fields wrote to 0/RANS" << endl;
  Info<< "cellDecompositionR list  wrote to RANS/constant \n" << endl;

  Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
      << "  ClockTime = " << runTime.elapsedClockTime() << " s"
      << nl << endl;
  
  Info << "End \n" << endl;
  
  return 0;
}


// ************************************************************************* //
