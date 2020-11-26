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

Namespace
    Foam::dynamicModels

Description

Class
    Foam::dynamicModel

Description
    Abstract base class for dynamic resolution models.

SourceFiles
    dynamicModel.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicModel_H
#define dynamicModel_H

#include "resolutionModel.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"
#include "wallDist.H"
#include "midPoint.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dynamicModel Declaration
\*---------------------------------------------------------------------------*/

class dynamicModel
:
    public resolutionModel,
    public IOdictionary
{

protected:

    // Protected data

        //- dict
        dictionary paramDict_;

        //- print parameters
        Switch printParams_;
        Switch printBounds_;

       // interpolated flag to find LES/RANS interface faces.
        midPoint< scalar > flagInterpolate_;
       // Distance to wall. To measure the distribution of the RANS cells
       wallDist distToWall_;

       //- resolution indicator. different meaning for different criteria
       volScalarField resIndicator_;

    // Protected member functions

        //- Print model coefficients
        virtual void printParams();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicModel(const dynamicModel&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicModel&);


public:

    //- Runtime type information
    TypeName("dynamicModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            dynamicModel,
            dictionary,
            (
                const fvMesh& meshL,
                const fvMesh& meshR
            ),
            (meshL, meshR)
        );


    // Constructors

        //- Construct from components
        dynamicModel
        (
            const word& type,
            const fvMesh& meshL,
            const fvMesh& meshR
        );


    // Selectors

        //- Return a reference to the selected RAS model
        static autoPtr<dynamicModel> New
        (
            const fvMesh& meshL,
            const fvMesh& meshR
        );


    //- Destructor
    virtual ~dynamicModel()
    {}


    // Member Functions

        // Access
        const dictionary& paramDict() const 
        {
            return paramDict();
        }
        

        //- Solve the turbulence equations and evaluate the turbulence viscosity
        virtual void evaluateFlag() = 0;

        //- Read DynamicProperties dictionary
        virtual bool read() = 0;

        //- Assertion
        //- Assert UpperLimit > LowerLimit
        void  assertOrder(scalar& l, scalar& u);
        void assertFlagInRange(volScalarField & lesFlagX);

        //- Diagnosis 
        virtual void printResolutionBounds() {};
        void diagnosisFlagConnectivity();
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
