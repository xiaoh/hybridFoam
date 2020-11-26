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
    Foam::incompressible::LESModels::muImpactADMRT

Description
    General base class for all incompressible models that can be implemented
    as an eddy viscosity, i.e. algebraic and one-equation models.

    Contains fields for k (SGS turbulent kinetic energy), gamma
    (modelled viscosity) and epsilon (SGS dissipation).

SourceFiles
    muImpactADMRT.C

\*---------------------------------------------------------------------------*/

#ifndef muImpactADMRT_H
#define muImpactADMRT_H

#include "muLESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class muImpactADMRT Declaration
\*---------------------------------------------------------------------------*/

class muImpactADMRT
:
 public muLESModel
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        muImpactADMRT(const muImpactADMRT&);
        muImpactADMRT& operator=(const muImpactADMRT&);

protected:

       volScalarField epsilonADM_;
       volScalarField k_;

public:

    TypeName("muImpactADMRT");

    // Constructors

        //- Construct from components
        muImpactADMRT
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const volSymmTensorField & RAvg,
            const volScalarField & epsilonAvg
        );


    //- Destructor
    virtual ~muImpactADMRT()
    {}


    // Member Functions

        //- Return SGS kinetic energy
       virtual tmp<volScalarField> k() const {return k_;}

        //- Return sub-grid disipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilonADM_;
        }

        //- Return the SGS viscosity.
        virtual tmp<volScalarField> nuSgs() const
        {
            return nu()*0.1;
        }

        virtual tmp<volScalarField> nuEff() const
        {
            FatalErrorIn("muImpactADMRT::nuEff()")
                << "This function is dummy and should never be called."
                << exit(FatalError);

            return tmp<volScalarField>
            (
                new volScalarField("nuEff", nuSgs() + nu())
            );
        }

        //- Return the sub-grid stress tensor.
       virtual tmp<volSymmTensorField> B() const { return ((2.0/3.0)*I)*k() - nuSgs()*twoSymm(fvc::grad(U())); }

        //- Return the effective sub-grid turbulence stress tensor
        //  including the laminar stress
        virtual tmp<volSymmTensorField> devBeff() const;

        //- Return the deviatoric part of the effective sub-grid
        //  turbulence stress tensor including the laminar stress
        virtual tmp<fvVectorMatrix> divDevBeff(volVectorField& U) const;

        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

       //- Enforce turbulence fields (relax forcing)
       virtual  void enforceFields
       (dimensionedScalar turbRelaxTime, scalar rampQ, word enforceMode) {};

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
