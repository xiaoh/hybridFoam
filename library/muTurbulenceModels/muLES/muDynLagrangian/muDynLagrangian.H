/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::incompressible::LESModels::muDynLagrangian

Description
    Lagrangian Two Equations Eddy Viscosity Model for incompressible
    flows

    \verbatim
        B = 2/3*k*I - 2*nuSgs*dev(D)
        Beff = 2/3*k*I - 2*nuEff*dev(D)

    where

        D = symm(grad(U))
        nuSgs = (flm/fmm)*delta^2*sqrt(2)*|D|
        nuEff = nuSgs + nu

        Two relaxation equations are used to evaluate flm and fmm:

        d/dt(flm) + div(U*flm)
        =
        (1/T)*(L && M - flm)

        d/dt(fmm) + div(U*fmm)
        =
        (1/T)*(M && M - flm)

    where

        L = F(U.U) - F(U).F(U)
        M = 2.0 delta^2 (F(|D|.D) - 4 F(|D|).F(D))
        T = 1.5*delta*(flm.fmm)^(-1.0/8.0)

    \endverbatim

    Reference:
    1. Charles Meneveau, Thomas Lund and William Cabot
    "A Lagrangian dynamic subgrid-scale model of turbulence"
    J. Fluid Mech (1996), vol 319, pp. 353-385

SourceFiles
    muDynLagrangian.C

\*---------------------------------------------------------------------------*/

#ifndef muDynLagrangian_H
#define muDynLagrangian_H

#include "muGenEddyVisc.H"
#include "simpleFilter.H"
#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                       Class muDynLagrangian Declaration
\*---------------------------------------------------------------------------*/

class muDynLagrangian
:
    public muGenEddyVisc
{
    // Private data

        volScalarField flm_;
        volScalarField fmm_;

        dimensionedScalar theta_;

        simpleFilter simpleFilter_;
        autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;

        dimensionedScalar flm0_;
        dimensionedScalar fmm0_;


    // Private Member Functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields
        (
            const tmp<volTensorField>& gradU
        );

        // Disallow default bitwise copy construct and assignment
        muDynLagrangian(const muDynLagrangian&);
        muDynLagrangian& operator=(const muDynLagrangian&);


public:

    //- Runtime type information
    TypeName("muDynLagrangian");

    // Constructors

        //- Construct from components
        muDynLagrangian
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const volSymmTensorField & RAvg,
            const volScalarField & epsilonAvg
        );


    //- Destructor
    virtual ~muDynLagrangian()
    {}


    // Member Functions

        //- Return SGS kinetic energy
        tmp<volScalarField> k(const tmp<volTensorField>& gradU) const
        {
            return 2.0*sqr(delta())*magSqr(dev(symm(gradU)));
        }

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k(fvc::grad(U()));
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nuSgs_ + nu())
            );
        }

       //- Enforce turbulence fields (relax forcing)
       virtual  void enforceFields
       (dimensionedScalar turbRelaxTime, scalar rampQ, word enforceMode) {};

        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);


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
