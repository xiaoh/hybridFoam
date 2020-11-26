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
    Foam::passiveModels::LKE

Description
    Wall distance model for turbulence resolution.

    The default model coefficients correspond to the following:
    @verbatim
        LKEParams
        {
    
        }
    @endverbatim

SourceFiles
    LKE.C

\*---------------------------------------------------------------------------*/

#ifndef LKE_H
#define LKE_H

#include "passiveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace passiveModels
{

/*---------------------------------------------------------------------------*\
                           Class LKE Declaration
\*---------------------------------------------------------------------------*/

class LKE
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
       volScalarField deltaL_;
       volScalarField kTotal_;
   
        //- Model coefficients
        const scalar volLengthWeight_;

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

        dimensionedScalar epsSmall_;

        //- compute delta
        void computeEquivDelta();

public:

    //- Runtime type information
    TypeName("LKE");

    // Constructors

        //- Construct from components
        LKE
        (
            const fvMesh& meshL,
            const fvMesh& meshR
        );


    //- Destructor
    virtual ~LKE()
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
