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

  Description
  A class for the interpolations between LES and RANS fields.

  SourceFiles
  meshTalkChannel.H
  meshTalkChannel.C
  \*---------------------------------------------------------------------------*/

#ifndef   meshTalkChannel_H
#define   meshTalkChannel_H

#include "meshToMesh.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  /*---------------------------------------------------------------------------* \
    Class meshTalkChannel Declaration
  \*---------------------------------------------------------------------------*/

  class meshTalkChannel
  {
    // Private data
  private:
      const fvMesh& meshL_;
      const fvMesh& meshR_;
      meshToMesh L2RInterp_;
      meshToMesh R2LInterp_;
      bool singleMesh_;
      
      const volScalarField & epsilonAvg_;
      const volVectorField & UAvg_;
      const  volSymmTensorField& RAvg_;
      const volScalarField & kR_;
      const volVectorField & UR_;
      const  volSymmTensorField& RR_;

      // defined on RANS mesh, mirrored from LES mesh
      volScalarField& mrEpsilonAvg_;
      volSymmTensorField& mrRAvg_;
      volVectorField mrUAvg_;

      // defined on LES mesh, mirrored from RANS mesh
      volScalarField mlkR_;
      volVectorField mlUR_;
      volSymmTensorField mlRR_;

      // indicator for cells which are out side of its partner
      labelIOList      outCellsL_;
      labelIOList      outCellsR_;

      // Private Member Functions

  public:

    // Constructors

    //- Construct from meshes
    meshTalkChannel( 
        const fvMesh& meshL,   const fvMesh& meshR ,
        volScalarField& mrEpsilonAvg,
        volSymmTensorField& mrRAvg
    );
    

    // Destructor
    ~meshTalkChannel();

    // Member Functions (public)
    //  template<template<class> class CombineOp>
    void interpolate(word mode="both");

    void smoothOverhangCellsL(); // as name suggest (via neighbour-averaging)
    void smoothOverhangCellsR();

    // Access
      const volVectorField & UAvg() {return singleMesh_? UAvg_ : mrUAvg_; }
      const volScalarField & epsilonAvg() { return singleMesh_? epsilonAvg_ : mrEpsilonAvg_;}
      const volSymmTensorField& RAvg()  { return singleMesh_? RAvg_ : mrRAvg_;}

    // defined on LES mesh, mirrored from RANS mesh
      const volVectorField & UR() { return singleMesh_? UR_ : mlUR_;}
      const volSymmTensorField & RR() {  return singleMesh_? RR_ : mlRR_;}
      const volScalarField & kR() {  return singleMesh_? kR_ : mlkR_;}

      
    
    // Check

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
