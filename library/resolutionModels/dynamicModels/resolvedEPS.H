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
    Foam::dynamicModels::resolvedEPS

Description
    Wall distance model for turbulence resolution.

    The default model coefficients correspond to the following:
    @verbatim
        resolvedEPSParams
        {
    
        }
    @endverbatim

SourceFiles
    resolvedEPS.C

\*---------------------------------------------------------------------------*/

#ifndef resolvedEPS_H
#define resolvedEPS_H

#include "dynamicModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dynamicModels
{

/*---------------------------------------------------------------------------*\
                           Class resolvedEPS Declaration
\*---------------------------------------------------------------------------*/

class resolvedEPS
:
    public dynamicModel
{
    // Protected data (to be inherited by Lagrangian averaging)
protected:
        //- Average (consistent) velocity
       const volVectorField& UAvg_;
       const volScalarField& epsilonAvg_;
       const volScalarField& nuSgs_;

       const IOdictionary relaxDict_;
       const IOdictionary transportDict_;

       const Time & runTime_;
       dimensionedScalar averagingTime_;
       dimensionedScalar nu_;

       //- Modeled (SGS) stress
       volScalarField epsilonMMAvg_;
   
        //- Model coefficients
        dimensionedScalar lower_;
        dimensionedScalar upper_;

        dimensionedScalar epsSmall_;

public:

    //- Runtime type information
    TypeName("resolvedEPS");

    // Constructors

        //- Construct from components
        resolvedEPS
        (
            const fvMesh& meshL,
            const fvMesh& meshR
        );


    //- Destructor
    virtual ~resolvedEPS()
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
