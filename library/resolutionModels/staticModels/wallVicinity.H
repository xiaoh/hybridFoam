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

Class
    Foam::staticModels::wallVicinity

Description
    Wall distance model for turbulence resolution.

    The default model coefficients correspond to the following:
    @verbatim
        wallVicinityParams
        {
    
        }
    @endverbatim

SourceFiles
    wallVicinity.C

\*---------------------------------------------------------------------------*/

#ifndef wallVicinity_H
#define wallVicinity_H

#include "staticModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace staticModels
{

/*---------------------------------------------------------------------------*\
                           Class wallVicinity Declaration
\*---------------------------------------------------------------------------*/

class wallVicinity
:
    public staticModel
{
    // Private data

        // Model coefficients
        dimensionedScalar ransRange_;

public:

    //- Runtime type information
    TypeName("wallVicinity");

    // Constructors

        //- Construct from components
        wallVicinity
        (
            const fvMesh& meshL,
            const fvMesh& meshR
        );


    //- Destructor
    virtual ~wallVicinity()
    {}


    // Member Functions

        //- evaluate resolution and update the LES/RANS flag
        virtual void evaluateFlag(){};

        void interpLesFlag();

        //- Read staticProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace staticModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
