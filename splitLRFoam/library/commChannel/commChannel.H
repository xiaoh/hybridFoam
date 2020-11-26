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
  commChannel

  Description
  A class for interpolation between LES and RANS in the context of hybrid method.
  Written from OF utility MapFields.

  SourceFiles
  commChannel.H
  commChannel.C
  \*---------------------------------------------------------------------------*/


#ifndef   commChannel_H
#define   commChannel_H

#include "IOdictionary.H"
#include "fvCFD.H"
#include "word.H"
#include "Time.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  /*---------------------------------------------------------------------------* \
    Class commChannel Declaration
  \*---------------------------------------------------------------------------*/

  class commChannel
  {
    // Private data

  private:
      Time& runTime_;
      fileName rootDirTarget_;
      fileName caseDirTarget_;

      fileName casePath_;
      fileName rootDirSource_;
      fileName caseDirSource_;
      word prefix_;
      HashTable<word> fieldsToMap_;
      HashSet<word> fieldsToReconstruct_;
      HashSet<word> fieldsToDecompose_;

      Time runTimeSource_;
      Time runTimeTarget_;

      HashTable<word> patchMap_;
      wordList cuttingPatches_;
 
    // Private Member Functions

      void mapConsistentMesh
      (
          const fvMesh& meshSource,
          const fvMesh& meshTarget
      );


  public:

    // Constructors
    
    //- Construct from components
    commChannel
    (
        Time& runTime,
        fileName& mappingSource,
        HashTable<word> fieldsToMap,
        HashSet<word> fieldsToReconstruct,
        HashSet<word> fieldsToDecompose
    );

    // Destructor
    ~commChannel(){};

    // Member Functions (public)
    void fetch(scalar sourceValue);
    void reconstruct();
    void decompose();

    // Access

    // Check
    static void printVersion();
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
