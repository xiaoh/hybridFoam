/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::homogeneousSpanIndex

Description
    Does averaging of fields over layers of cells. Assumes layered mesh.

SourceFiles
    homogeneousSpanIndex.C

\*---------------------------------------------------------------------------*/

#ifndef homogeneousSpanIndex_H
#define homogeneousSpanIndex_H

#include "regionSplit.H"
#include "direction.H"
#include "scalarField.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class homogeneousSpanIndex Declaration
\*---------------------------------------------------------------------------*/

class homogeneousSpanIndex
{

    // Private data
  
        const vector spanDirection_;

        //- Per cell the global region
        autoPtr<regionSplit> cellRegion_;

        //- Per global region the number of cells (scalarField so we can use
        //  field algebra)
        scalarField regionCount_;


    // Private Member Functions

        void calcLayeredRegions
        (
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construct and assignment
        homogeneousSpanIndex(const homogeneousSpanIndex&);
        void operator=(const homogeneousSpanIndex&);


public:

    // Constructors

        //- Construct from dictionary
        homogeneousSpanIndex(const fvMesh&);

    // Member Functions

        // Access

            //- Sum field per region
            template<class T>
            Field<T> regionSum(const Field<T>& cellField) const;

            //- collapse a field to a line
            template<class T>
            Field<T> collapse
            (
                const Field<T>& vsf
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "homogeneousSpanIndexTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
