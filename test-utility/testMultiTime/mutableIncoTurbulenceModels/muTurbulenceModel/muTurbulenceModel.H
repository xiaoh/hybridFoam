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
    Foam::incompressible::muTurbulenceModels

Description
    Namespace for incompressible turbulence turbulence models.

Class
    Foam::incompressible::muTurbulenceModel

Description
    Abstract base class for incompressible turbulence models
    (RAS, LES and laminar).

SourceFiles
    muTurbulenceModel.C
    newTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef muTurbulenceModel_H
#define muTurbulenceModel_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "incompressible/transportModel/transportModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;

namespace incompressible
{

/*---------------------------------------------------------------------------*\
                           Class muTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

class muTurbulenceModel
{

protected:

    // Protected data

        const Time& runTime_;
        const fvMesh& mesh_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        transportModel& transportModel_;

  //- Reference to the stress/dissipation/flag fields in relaxForcing
        const volSymmTensorField & RAvg_;
        const volScalarField & epsilonAvg_;
        const volScalarField & lesFlag_;



private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        muTurbulenceModel(const muTurbulenceModel&);

        //- Disallow default bitwise assignment
        void operator=(const muTurbulenceModel&);


public:

    //- Runtime type information
    TypeName("muTurbulenceModel");


    // Declare run-time New selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            muTurbulenceModel,
            muTurbulenceModel,
            (
                const volVectorField& U,
                const surfaceScalarField& phi,
                transportModel& lamTransportModel,
                const volSymmTensorField & RAvg,
                const volScalarField & epsilonAvg,
                const volScalarField & lesFlag
            ),
            (U, phi, lamTransportModel, RAvg, epsilonAvg, lesFlag)
        );


    // Constructors

        //- Construct from components
        muTurbulenceModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& lamTransportModel,
            const volSymmTensorField & RAvg,
            const volScalarField & epsilonAvg,
            const volScalarField & lesFlag
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<muTurbulenceModel> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& lamTransportModel,
            const volSymmTensorField & RAvg,
            const volScalarField & epsilonAvg,
            const volScalarField & lesFlag
        );


    //- Destructor
    virtual ~muTurbulenceModel()
    {}


    // Member Functions

        //- Access function to velocity field
        inline const volVectorField& U() const
        {
            return U_;
        }

        //- Access function to flux field
        inline const surfaceScalarField& phi() const
        {
            return phi_;
        }

        //- Access function to incompressible transport model
        inline transportModel& transport() const
        {
            return transportModel_;
        }

        //- Return the laminar viscosity
        const volScalarField& nu() const
        {
            return transportModel_.nu();
        }

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const = 0;

        //- Return the effective viscosity
        virtual tmp<volScalarField> nuEff() const = 0;

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const = 0;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const = 0;

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const = 0;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const = 0;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const = 0;

       //- Enforce turbulence fields (relax forcing)
        virtual  void enforceFields
         (dimensionedScalar turbRelaxTime, scalar rampQ, word enforceMode) = 0;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct() = 0;

        //- Read turbulenceProperties dictionary
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
