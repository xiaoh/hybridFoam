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
    Foam::dynamicModels::resolvedTKE

Description
    Wall distance model for turbulence resolution.

    The default model coefficients correspond to the following:
    @verbatim
        resolvedTKEParams
        {
    
        }
    @endverbatim

SourceFiles
    resolvedTKE.C

\*---------------------------------------------------------------------------*/

#ifndef resolvedTKE_H
#define resolvedTKE_H

#include "dynamicModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dynamicModels
{

/*---------------------------------------------------------------------------*\
                           Class resolvedTKE Declaration
\*---------------------------------------------------------------------------*/

class resolvedTKE
:
    public dynamicModel
{
    // Protected data (to be inherited by Lagrangian averaging)
protected:
        //- Resolved stress
        const volSymmTensorField& RResolvedAvg_;

        //- Modeled (SGS) stress
        const volScalarField& kSgsAvg_;

        //- turb/cell length scale fields
       volScalarField kTotal_;
   
        //- Model coefficients
        dimensionedScalar lower_;
        dimensionedScalar upper_;

        dimensionedScalar epsSmall_;

public:

    //- Runtime type information
    TypeName("resolvedTKE");

    // Constructors

        //- Construct from components
        resolvedTKE
        (
            const fvMesh& meshL,
            const fvMesh& meshR
        );


    //- Destructor
    virtual ~resolvedTKE()
    {}


    // Member Functions

        //- print bounds
        virtual void printResolutionBounds();

        //- evaluate resolution and update the LES/RANS flag
        virtual void evaluateFlag();

        //- Read dynamicProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dynamicModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
