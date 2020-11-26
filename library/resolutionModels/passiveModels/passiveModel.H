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
    Foam::passiveModels

Description

Class
    Foam::passiveModel

Description
    Abstract base class for passive resolution models.
    These models are used for the investigation of resolution criteria, 
    so that several criteria can be used in one single run.
    They do not set lesFlagL and lesflagR, and should NOT be used in actual hybrid simulations.

SourceFiles
    passiveModel.C

\*---------------------------------------------------------------------------*/

#ifndef passiveModel_H
#define passiveModel_H

#include "resolutionModel.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class passiveModel Declaration
\*---------------------------------------------------------------------------*/

class passiveModel
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

    // Protected member functions

        //- Print model coefficients
        virtual void printParams();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        passiveModel(const passiveModel&);

        //- Disallow default bitwise assignment
        void operator=(const passiveModel&);


public:

    //- Runtime type information
    TypeName("passiveModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            passiveModel,
            dictionary,
            (
                const fvMesh& meshL,
                const fvMesh& meshR
            ),
            (meshL, meshR)
        );


    // Constructors

        //- Construct from components
        passiveModel
        (
            const word& type,
            const fvMesh& meshL,
            const fvMesh& meshR
        );


    // Selectors

        //- Return a reference to the selected RAS model
        static autoPtr<passiveModel> New
        (
            const fvMesh& meshL,
            const fvMesh& meshR
        );


    //- Destructor
    virtual ~passiveModel()
    {}


    // Member Functions

        // Access
        const dictionary& paramDict() const 
        {
            return paramDict();
        }
        

        //- Solve the turbulence equations and evaluate the turbulence viscosity
       virtual void evaluateFlag() = 0;

        //- Read PassiveProperties dictionary
        virtual bool read() = 0;

       //- Diagnosis 
        virtual void printResolutionBounds() {};
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
