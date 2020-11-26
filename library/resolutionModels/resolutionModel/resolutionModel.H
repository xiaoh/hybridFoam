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
    Foam::resolutionModels

Description
    Namespace for turbulence resolution models.

Class
    Foam::resolutionModel

Description
    Abstract base class for resolution models
    (static: wall distance; 
    dynamic: energy, length, dissipation).

SourceFiles
    resolutionModel.C

\*---------------------------------------------------------------------------*/

#ifndef resolutionModel_H
#define resolutionModel_H

#include "volFieldsFwd.H"
#include "volFields.H"
#include "Time.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "meshToMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class resolutionModel Declaration
\*---------------------------------------------------------------------------*/

class resolutionModel
{

protected:

    // Protected data
    
        const Time& runTime_;
        const fvMesh& meshL_;
        const fvMesh& meshR_;
    
        label nLocal_; // no. of cells in local LES mesh
        label nGlobal_; // no. of cells in global LES mesh
    
    // This is what the whole class is about
        volScalarField lesFlagL_; // Interpolated to LES mesh
        volScalarField lesFlagR_; // Computed on RANS mesh

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        resolutionModel(const resolutionModel&);

        //- Disallow default bitwise assignment
        void operator=(const resolutionModel&);


public:

    //- Runtime type information
    TypeName("resolutionModel");


    // Declare run-time New selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            resolutionModel,
            resolutionModel,
            (
                const fvMesh& meshL,
                const fvMesh& meshR
            ),
            (meshL, meshR)
        );


    // Constructors

        //- Construct from components
        resolutionModel
        (
            const fvMesh& meshL,
            const fvMesh& meshR
        );


    // Selectors

        //- Return a reference to the selected resolution model
        static autoPtr<resolutionModel> New
        (
            const fvMesh& meshL,
            const fvMesh& meshR
        );


    //- Destructor
    virtual ~resolutionModel()
    {}


    // Member Functions

        //- return lesFlag on LES mesh
        virtual const volScalarField& lesFlagL() const
          {return lesFlagL_;}

        //- return lesFlag on RANS mesh
        virtual const volScalarField& lesFlagR() const
          {return lesFlagR_;}

        void interpLesFlagL2R();
        void interpLesFlagR2L();

        //- Diagnosis LES region
        void diagnosisFlagRegion();

        //- Evaluate resolution
        virtual void evaluateFlag() = 0;

        //- Read turbulenceProperties dictionary
        virtual bool read() = 0;

        static void printVersion();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
