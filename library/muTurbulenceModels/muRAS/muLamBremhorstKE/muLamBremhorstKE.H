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
    Foam::incompressible::RASModels::muLamBremhorstKE

Description
    Lam and Bremhorst low-Reynolds number k-epsilon turbulence model
    for incompressible flows

SourceFiles
    muLamBremhorstKE.C

\*---------------------------------------------------------------------------*/

#ifndef muLamBremhorstKE_H
#define muLamBremhorstKE_H

#include "muRASModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class muLamBremhorstKE Declaration
\*---------------------------------------------------------------------------*/

class muLamBremhorstKE
:
    public muRASModel
{
    // Private data

        dimensionedScalar Cmu_;
        dimensionedScalar C1_;
        dimensionedScalar C2_;
        dimensionedScalar sigmaEps_;

        volScalarField k_;
        volScalarField epsilon_;

        wallDist y_;
        volScalarField Rt_;

        volScalarField fMu_;
        volScalarField nut_;

        // Forcing fields
        volScalarField Qk_;
        volScalarField Qepsilon_;
        volScalarField kLES_;

        void diagnosisDiffKE(string msg) const;

public:

    //- Runtime type information
    TypeName("muLamBremhorstKE");


    // Constructors

        //- Construct from components
        muLamBremhorstKE
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const volSymmTensorField & RAvg,
            const volScalarField & epsilonAvg
        );


    //- Destructor
    virtual ~muLamBremhorstKE()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nut_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_/sigmaEps_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

       //- Enforce turbulence fields via relax forcing
        virtual  void enforceFields
        (dimensionedScalar turbRelaxTime, scalar rampQ, word enforceMode);

        virtual void diagnosis() const;

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
