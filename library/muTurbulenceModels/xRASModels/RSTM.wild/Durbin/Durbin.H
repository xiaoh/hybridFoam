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
    Foam::incompressible::RASModels::Durbin

Description
    Durbin Elliptic Relaxation Model.

    Durbin, A Reynolds stress model for near-wall turbulence, J. Fluid
    Mech. vol 249, 1993.

    with the possibility to turn on the ellliptic operator correction
    by Manceau et. al. (2001):

    Manceau, Wang, and Laurence. Inhomogeneity and anisotropy effects
    on the redistribution term in Reynolds-averaged Navier-Stokes
    modelling. J. Fluid Mech. vol 438, 2001.

    The default model coefficients correspond to the following:
    @verbatim
        DurbinCoeffs
        {
            yStarLim        5;
            CL              0.2;
            CmuKE           0.09;
            Clrr1           1.22;
            Clrr2           0.6;
            C1              1.44;
            C2              1.9;
            Cmu             0.23;
            sigmaEps        1.65;
            sigmaK          1.2;
            CEta            80;
            solveK          true;
            SSG             false;
            Cg1             3.4;
            Cg1s            1.8;
            Cg2             4.2;
            Cg3             0.8;
            Cg3s            1.3;
            Cg4             1.25;
            Cg5             0.4;
            ellipticOperatorCorrection false;
            beta            0.083;
         }
    @endverbatim

SourceFiles
    Durbin.C

\*---------------------------------------------------------------------------*/

#ifndef Durbin_H
#define Durbin_H

#include "GenElliptic.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                            Class Durbin Declaration
\*---------------------------------------------------------------------------*/

class Durbin
:
    public GenElliptic
{
    // Private data

        // coeffs related to elliptic relaxation
        Switch solveK_;

        Switch ellipticOperatorCorrection_;

        //Manceau correction related coefficient
        dimensionedScalar beta_;

        // Fields

            volSymmTensorField f_;
public:

    //- Runtime type information
    TypeName("Durbin");

    // Constructors

        //- Construct from components
        Durbin
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    //- Destructor
    virtual ~Durbin()
    {}


    // Member Functions

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

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
