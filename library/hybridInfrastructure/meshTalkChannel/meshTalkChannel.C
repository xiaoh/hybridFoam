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
  meshTalkChannel

  \*---------------------------------------------------------------------------*/

#include "meshTalkChannel.H"

namespace Foam {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
meshTalkChannel::meshTalkChannel
(
    const fvMesh& meshL, const fvMesh& meshR, 
    volScalarField& mrEpsilonAvg,
    volSymmTensorField& mrRAvg
)
    :
    meshL_(meshL), // this is a reference to a fvMesh object in *.H: const fvMesh& meshL_;
    meshR_(meshR),
    L2RInterp_(meshL, meshR, HashTable<word>(0), wordList()),// this is a meshToMesh object, defined in *.H
    // meshToMesh is defined in meshToMesh.H in the OF2.2.0 original version, constructor declaration:
        // meshToMesh
        // (
        //     const fvMesh& fromMesh,
        //     const fvMesh& toMesh,
        //     const HashTable<word>& patchMap,
        //     const wordList& cuttingPatchNames
        // );
    R2LInterp_(meshR, meshL, HashTable<word>(0), wordList()),// this is a meshToMesh object, defined in *.H
    singleMesh_(meshL==meshR),

    epsilonAvg_( meshL.lookupObject<const volScalarField>("epsilonAvg") ),
    UAvg_( meshL.lookupObject<const volVectorField>("UAvg") ),
    RAvg_( meshL.lookupObject<const volSymmTensorField>("RAvg") ),
    
    kR_( meshR.lookupObject<const volScalarField>("kR") ),
    UR_( meshR.lookupObject<const volVectorField>("UR") ),
    RR_( meshR.lookupObject<const volSymmTensorField>("RR") ),

    // Fields defined on LES mesh but mirrored to RANS mesh
    // LES --> RANS

    mrEpsilonAvg_(mrEpsilonAvg),
    mrRAvg_(mrRAvg),

    mrUAvg_
      (
          IOobject
          (
              "mrUAvg",
              meshR.time().timeName(),
              meshR,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE
          ),
          UR_
      ),

    // Fields defined on RANS mesh but mirrored to LES mesh
    // RANS --> LES 

    mlkR_
    (
     IOobject
     (
      "mlkR",
      meshR.time().timeName(),
      meshL,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
      ),
     meshL,
     dimensionedScalar
     ("mlkR", dimVelocity*dimVelocity, scalar(0.0))
    ),

    mlUR_
      (
          IOobject
          (
              "mlUR",
              meshR.time().timeName(),
              meshL,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE
          ),
          meshL, 
          dimensionedVector
          ("mlUR", dimVelocity, vector::zero)
      ),
    // Reynolds stress from RANS model
    mlRR_
      (
          IOobject
          (
              "mlRR",
              meshL.time().timeName(),
              meshL,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE
          ),
          meshL,
          dimensionedSymmTensor
          ("mlRR", dimVelocity*dimVelocity, symmTensor::zero)
      ),

    outCellsL_
    (
        IOobject
        (
            "overHangCellsL",
            meshL.time().constant(),
            meshL,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        0
    ),
    
    outCellsR_
    (
        IOobject
        (
            "overHangCellsR",
            meshR.time().constant(),
            meshR,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        0
    )
    
{   
    if(singleMesh_)
    {
        Info << "NOTE: Identical LES and RANS meshes. No interpolations." << endl;
        mrUAvg_.writeOpt() = IOobject::NO_WRITE;
        mlUR_.writeOpt() = IOobject::NO_WRITE;
        mlkR_.writeOpt() = IOobject::NO_WRITE;
        mlRR_.writeOpt() = IOobject::NO_WRITE;
        return;
    }
    
    if(!Pstream::parRun()) // serial run
    {
        Info << "Serial run. Skipping overhang cell list." << endl;
    }
    else
    {
        if (! returnReduce(outCellsL_.headerOk(), andOp<bool>()) )
        {
            DynamicList<label> tmpL(256);
            const  volVectorField & cellCenters = meshL.C();
            Info << "LES: Over-hang cell list not found. computing ..." << endl;
            forAll(cellCenters, celli)
            {
                point p = cellCenters[celli];
                if(meshR.findCell(p) < 0) // Find L cell in RANS domain
                {
                    tmpL.append(celli); 
                }
            }
            outCellsL_.transfer(tmpL);
            outCellsL_.write();
        }
        else
        {
            Info << "LES: Over-hang cell list read successfully." << endl;
        }
        
        if (! returnReduce(outCellsR_.headerOk(), andOp<bool>() ))
        {
            DynamicList<label> tmpR(256);
            const  volVectorField & cellCenters = meshR.C();
            Info << "RAS: Over-hang cell list not found. computing ..." << endl;
            forAll(cellCenters, celli)
            {
                point p = cellCenters[celli];
                if(meshL.findCell(p) < 0) // Find R cell in LES domain
                {
                    tmpR.append(celli); 
                }
            }
            outCellsR_.transfer(tmpR);
            outCellsR_.write();
        }
        else
        {
            Info << "RAS: Over-hang cell list read successfully." << endl;
        }
    }

    interpolate("both");
}
  
  // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

  meshTalkChannel::~meshTalkChannel()
  {}

  // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
  
  //- Interpolate both ways (LES->RANS, RANS->LES)


  //  template<template<class> class CombineOp>
  void meshTalkChannel::interpolate(word mode) 
  {
      if(singleMesh_) return;

      bool L2Rflag = false;
      bool R2Lflag = false;

      if(mode == "L2R")
      {
          L2Rflag = true;
      }
      else if(mode == "R2L")
      {
          R2Lflag = true;
      }
      else
      {
          R2Lflag = true;
          L2Rflag = true;
      }

      //- LES --> RANS
      if(L2Rflag)
      {
          L2RInterp_.interpolateInternalField
              (
	       mrUAvg_, UAvg_, 
	       meshToMesh::INTERPOLATE,
	       eqOp<vector>()
              );
          
          L2RInterp_.interpolateInternalField
              (
	       mrEpsilonAvg_, epsilonAvg_,  
	       meshToMesh::INTERPOLATE,
	       eqOp<scalar>()
              );
          
          L2RInterp_.interpolateInternalField
              (
	       mrRAvg_, RAvg_, 
	       meshToMesh::INTERPOLATE,
	       eqOp<symmTensor>()
              );
   
          if(outCellsR_.size() && Pstream::parRun())
          {
              smoothOverhangCellsR();
          }
      }
      
      //- RANS --> LES
      
      if(R2Lflag)
      {
          R2LInterp_.interpolateInternalField
              (
	       mlUR_, UR_, 
	       meshToMesh::INTERPOLATE,
	       eqOp<vector>()
              );
          
          R2LInterp_.interpolateInternalField
              (
	       mlRR_, RR_,
	       meshToMesh::INTERPOLATE,
	       eqOp<symmTensor>()
              );
          
          R2LInterp_.interpolateInternalField
              (
	       mlkR_, kR_,
	       meshToMesh::INTERPOLATE,
	       eqOp<scalar>()
              );
          
          if(outCellsL_.size() && Pstream::parRun())
          {
              smoothOverhangCellsL();
          }
      }
      
  }


void meshTalkChannel::smoothOverhangCellsL()
{
    labelList outstanding(outCellsL_); // copy
    DynamicList<label> leftover(0);
    labelList illCells(meshL_.nCells(), 0); // indicator field
    
    forAll(outstanding, ii) 
    {
        label celli = outstanding[ii];
        illCells[celli] = 1;
    }

    label nSmooth = 0;
    label nMax = 30;
    const labelListList& cellCell = meshL_.cellCells();
    while(outstanding.size())
    {
        forAll(outstanding, ii)
        {
            label celli = outstanding[ii];
            labelList neighbours = cellCell[celli];
            
            label nHealthyNei = 0;
            vector sumUR = vector::zero;
            double sumkR = 0.0;
            symmTensor sumRR = symmTensor::zero;
            for(label jj=0; jj < neighbours.size(); jj++)
            {
                label neiID = neighbours[jj];
                if(illCells[neiID] == 0) // neighbor is a healthy cell
                {
                    sumUR += mlUR_[neiID];
                    sumkR += mlkR_[neiID];
                    sumRR += mlRR_[neiID];
                    nHealthyNei++;
                }
            }

            if (! nHealthyNei) // all neighbours are ill: wait for next round of smoothing
            {
                leftover.append(celli);
            }
            else // celli is cured: set indicator field.
            {
                mlUR_[celli] = sumUR / nHealthyNei;
                mlkR_[celli] = sumkR / nHealthyNei;
                mlRR_[celli] = sumRR / nHealthyNei; 
                illCells[celli] = 0;
            }
        }

        nSmooth++;

        if(nSmooth >= nMax)
        {
            Pout << "Warning in meshTalkChannel::smoothOverhangCellsR: Cells outstanding after max attemps: " << nMax << endl;
            break;
        }
        else
        {
            outstanding = leftover;
            leftover.clearStorage();
        }     
    }
}

    
void meshTalkChannel::smoothOverhangCellsR()
{
    labelList outstanding(outCellsR_); // copy
    DynamicList<label> leftover(0);
    labelList illCells(meshR_.nCells(), 0); // indicator field
    
    forAll(outstanding, ii) 
    {
        label celli = outstanding[ii];
        illCells[celli] = 1;
    }

    label nSmooth = 0;
    label nMax = 30;
    const labelListList& cellCell = meshR_.cellCells();
    while(outstanding.size())
    {
        forAll(outstanding, ii)
        {
            label celli = outstanding[ii];
            labelList neighbours = cellCell[celli];
            
            label nHealthyNei = 0;
            vector sumUAvg = vector::zero;
            double sumEps = 0.0;
            symmTensor sumRAvg = symmTensor::zero;
            for(label jj=0; jj < neighbours.size(); jj++)
            {
                label neiID = neighbours[jj];
                if(illCells[neiID] == 0) // neighbor is a healthy cell
                {
                    sumUAvg += mrUAvg_[neiID];
                    sumEps  += mrEpsilonAvg_[neiID];
                    sumRAvg += mrRAvg_[neiID];
                    nHealthyNei++;
                }
            }

            if (! nHealthyNei) // all neighbours are ill: wait for next round of smoothing
            {
                leftover.append(celli);
            }
            else // celli is cured: set indicator field.
            {
                mrUAvg_[celli]        = sumUAvg / nHealthyNei;
                mrEpsilonAvg_[celli]  = sumEps  / nHealthyNei;
                mrRAvg_[celli]        = sumRAvg / nHealthyNei; 
                illCells[celli] = 0;
            }
        }

        nSmooth++;

        if(nSmooth >= nMax)
        {
            Pout << "Warning in meshTalkChannel::smoothOverhangCellsL: Cells outstanding after max attemps: " << nMax << endl;
            break;
        }
        else
        {
            outstanding = leftover;
            leftover.clearStorage();
        }     
    }

}



  // * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

  // * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

  // * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace Foam


// ************************************************************************* //
