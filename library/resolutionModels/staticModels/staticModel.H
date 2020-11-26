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
    Foam::staticModels

Description

Class
    Foam::staticModel

Description
    Abstract base class for static resolution models.

SourceFiles
    staticModel.C

\*---------------------------------------------------------------------------*/

#ifndef staticModel_H
#define staticModel_H

#include "resolutionModel.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class staticModel Declaration
\*---------------------------------------------------------------------------*/

class staticModel
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
        staticModel(const staticModel&);

        //- Disallow default bitwise assignment
        void operator=(const staticModel&);


public:

    //- Runtime type information
    TypeName("staticModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            staticModel,
            dictionary,
            (
                const fvMesh& meshL,
                const fvMesh& meshR
            ),
            (meshL, meshR)
        );


    // Constructors

        //- Construct from components
        staticModel
        (
            const word& type,
            const fvMesh& meshL,
            const fvMesh& meshR
        );


    // Selectors

        //- Return a reference to the selected RAS model
        static autoPtr<staticModel> New
        (
            const fvMesh& meshL,
            const fvMesh& meshR
        );


    //- Destructor
    virtual ~staticModel()
    {}


    // Member Functions

        // Access
        const dictionary& paramDict() const 
        {
            return paramDict();
        }
        

        //- Solve the turbulence equations and evaluate the turbulence viscosity
        virtual void evaluateFlag() = 0;

        //- Read StaticProperties dictionary
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
