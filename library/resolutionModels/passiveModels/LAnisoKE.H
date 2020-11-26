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
    Foam::passiveModels::LAnisoKE

Description
    Wall distance model for turbulence resolution.

    The default model coefficients correspond to the following:
    @verbatim
        LAnisoKEParams
        {
    
        }
    @endverbatim

SourceFiles
    LAnisoKE.C

\*---------------------------------------------------------------------------*/

#ifndef LANISOKE_H
#define LANISOKE_H

#include "passiveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace passiveModels
{

/*---------------------------------------------------------------------------*\
                           Class LAnisoKE Declaration
\*---------------------------------------------------------------------------*/

class LAnisoKE
:
    public passiveModel
{
    // Protected data (to be inherited by Lagrangian averaging)
protected:
        const Time & runTime_;

        //- Resolved stress
        const volSymmTensorField& RResolvedAvg_;

        //- Modeled (SGS) stress
        const volScalarField& kSgsAvg_;
    
       // average velocity (mean motion)
       const volVectorField& UAvg_;

        //- avg dissipation
        const volScalarField& epsilonAvg_;

        //- turb/cell length scale fields
       volScalarField turbLengthScale_;
       volScalarField magDeltaL_;
       volScalarField kTotal_;

       volVectorField turbLengthScaleV_;
       volVectorField deltaL_;
       volVectorField kTotalV_;
   
       //- Model coefficients
       const dimensionedVector scaleTurbLengthBy_;

       const volScalarField& nuSgs_;

       const IOdictionary relaxDict_;
       const IOdictionary transportDict_;

       dimensionedScalar averagingTime_;
       dimensionedScalar nu_;

       //- Modeled (SGS) stress
       volScalarField epsilonMMAvg_;

       //- resolution indicators for L, K, E
       volScalarField resIndicatorL_;
       volScalarField resIndicatorK_;
       volScalarField resIndicatorE_;

       //- resolution indicators for LAniso
       volScalarField resIndX_;
       volScalarField resIndY_;
       volScalarField resIndZ_;
       volScalarField resIndMin_;

        dimensionedScalar epsSmall_;

        //- compute delta
        void computeEquivDelta();

public:

    //- Runtime type information
    TypeName("LAnisoKE");

    // Constructors

        //- Construct from components
        LAnisoKE
        (
            const fvMesh& meshL,
            const fvMesh& meshR
        );


    //- Destructor
    virtual ~LAnisoKE()
    {}


    // Member Functions

        //- evaluate resolution and update the LES/RANS flag
        virtual void evaluateFlag();

        //- Read passiveProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace passiveModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
