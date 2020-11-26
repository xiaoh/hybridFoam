/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    muV2F

Description
    mu V2F turbulence model.

SourceFiles
    muV2F.C

\*---------------------------------------------------------------------------*/

#ifndef muV2F_H
#define muV2F_H

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
                           Class muV2F Declaration
\*---------------------------------------------------------------------------*/

class muV2F
:
    public muRASModel
{
    // Private data

        dimensionedScalar Cmu;
        dimensionedScalar CmuKE;
        dimensionedScalar Ceps10;
        dimensionedScalar Ceps11;
        dimensionedScalar Ceps2;
        dimensionedScalar C1;
        dimensionedScalar C2;
        dimensionedScalar CL;
        dimensionedScalar CEta;
        dimensionedScalar oneOnSigmaK;
        dimensionedScalar oneOnSigmaEps;
        dimensionedScalar yStarLim;

        dimensionedScalar f0_;

        wallDist yw_;

        volScalarField k_;
        volScalarField epsilon_;
        volScalarField v2_;
        volScalarField f_;

        volScalarField nut_;

       // Forcing fields
            volScalarField Qk_;
            volScalarField Qepsilon_;
            volScalarField kLES_;

    // Private member functions

        tmp<volScalarField> T() const;
        tmp<volScalarField> L() const;
        void diagnosisDiffKE(string msg) const;
        void checkEnforceMode(word enforceMode);
public:

    //- Runtime type information
    TypeName("muV2F");

    // Constructors

        //- from components
        muV2F
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const volSymmTensorField & RAvg,
            const volScalarField & epsilonAvg
        );


    // Destructor

        ~muV2F(){}


    // Member Functions

        tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField> 
            (
                new volScalarField("DkEff", oneOnSigmaK*nut_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField> 
            (
                new volScalarField("DepsilonEff", oneOnSigmaEps*nut_ + nu())
            );
        }

        tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Note that epsilonTilda is returned as epsilon.
        //  This is the appropriate variable for most purposes.
        tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        tmp<volScalarField> v2() const
        {
            return v2_;
        }

        tmp<volScalarField> f() const
        {
            return f_;
        }

        tmp<volSymmTensorField> R() const;

//         tmp<fvVectorMatrix> divR(volVectorField& U) const;
    
        //- Return the effective stress tensor including the laminar stress
        tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        void correct();

       //- Enforce turbulence fields via relax forcing
        virtual  void enforceFields
        (dimensionedScalar turbRelaxTime, scalar rampQ, word enforceMode);

        virtual void diagnosis() const;

        //- Read turbulenceProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
