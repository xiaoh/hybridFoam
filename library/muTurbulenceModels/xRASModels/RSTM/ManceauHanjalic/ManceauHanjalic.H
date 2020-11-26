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
    Foam::incompressible::RASModels::ManceauHanjalic

Description

    Elliptic blending model based on Manceau and Hanjalic (2002).
    Implementation: Heng Xiao (xiaoheng2011@gmail.com)
    Reivew: Michael Wild
    Bug report and comments to: HX.

    Part of the code is borrowed from the LRR model implentation in OF
    distribution and from the V2F model implementation.
    
    The basic idea of the model is according to:

    Manceau and Hanjalic, Elliptic blending model: a new
    near-wall Reynolds-stress turbulence closure. Physics of Fluids 14
    (74), 2002.
    
    However, the implementation follows the modified version in
   
    Thielen et al., Prediction of flow and heat transfer in multiple
    impinging jets with an elliptic-blending second moment closure.
    Int. J. Heat and Mass Transfer, vol 48 (8), 2004.

    The default model coefficients correspond to the following:
    @verbatim
        ManceauHanjalicCoeffs
    {
        CmuKE           0.09;
        Clrr1           1.22;
        Clrr2           0.6; 
        C1              1.44;
        C2              1.85;
        Cmu             0.22;
        sigmaEps        1.65; 
        CL              0.25; // one should use 0.161 instead for ManceauHanjalic model.
        CEta            80; 
        yStarLim        5;   
        SSG             false;
        Cg1             3.4;  // Cg1-5 effective if SSG is turned on.
        Cg1s            1.8;
        Cg2             4.2;
        Cg3             0.8;
        Cg3s            1.3;
        Cg4             1.25;
        Cg5             0.4; 
    }
    @endverbatim

SourceFiles
    ManceauHanjalic.C

\*---------------------------------------------------------------------------*/

#ifndef ManceauHanjalic_H
#define ManceauHanjalic_H

#include "GenElliptic.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class ManceauHanjalic Declaration
\*---------------------------------------------------------------------------*/

class ManceauHanjalic
:
    public GenElliptic
{
    // Private data
      // Fields
            volScalarField alpha_;
            Switch crossTurbDiffusion_; // turn on cross diffusion due to Reynolds stresses
public:

    //- Runtime type information
    TypeName("ManceauHanjalic");

    // Constructors

        //- Construct from components
        ManceauHanjalic
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    //- Destructor
    virtual ~ManceauHanjalic()
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
