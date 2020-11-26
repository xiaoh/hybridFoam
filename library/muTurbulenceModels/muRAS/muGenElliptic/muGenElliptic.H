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
    Foam::incompressible::RASModels::muGenElliptic

Description
    An abstract virtual class for elliptic relaxation RSTM models.
    Examples include:

    Durbin Elliptic Relaxation Model.
    
    Durbin, A Reynolds stress model for near-wall turbulence, J. Fluid
    Mech. vol 249, 1993.

    and 
    Hanjalic Elliptic Blending Model

    Manceau and Hanjalic, Elliptic blending model: a new
    near-wall Reynolds-stress turbulence closure. Physics of Fluids 14
    (74), 2002.
    
    The default model coefficients correspond to the following:
    @verbatim
        Coeffs
        {
          yStarLim        5;               
          CL              0.25;             
          CmuKE           0.09;            
          Clrr1           1.0;            
          Clrr2           0.6;             
          C1              1.44;            
          C2              1.9;             
          Cmu             0.23;            
          sigmaEps        1.65;            
          sigmaK          1.2;             
          CEta            80;              
          SSG             false;           
          Cg1             3.4;             
          Cg1s            1.8;             
          Cg2             4.2;             
          Cg3             0.8;             
          Cg3s            1.3;             
          Cg4             1.25;            
          Cg5             0.4;             
         }
    @endverbatim

SourceFiles
    muGenElliptic.C

\*---------------------------------------------------------------------------*/

#ifndef muGenElliptic_H
#define muGenElliptic_H

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
                           Class muGenElliptic Declaration
\*---------------------------------------------------------------------------*/

class muGenElliptic
:
    virtual public muRASModel
{
    // Private member functions:

    // Disallow default bitwise copy construct and assignment
    muGenElliptic(const muGenElliptic&); 
    muGenElliptic& operator=(const muGenElliptic&); 

    //Protected data

    protected:
        // Model coefficients
            const fvMesh& mesh_;
            dimensionedScalar CmuKE_;

            // LRR-IP model related coeffs
            dimensionedScalar Clrr1_;
            dimensionedScalar Clrr2_;

            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar Cmu_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar sigmaK_;

            // coeffs related to elliptic relaxation
            dimensionedScalar CL_;
            dimensionedScalar CEta_;
            dimensionedScalar yStarLim_; 
            const Switch      implicitDiv_;

            // SSG model related coeffs
            const Switch  SSG_;
            dimensionedScalar Cg1_;
            dimensionedScalar Cg1s_;
            dimensionedScalar Cg2_;
            dimensionedScalar Cg3_;
            dimensionedScalar Cg3s_;
            dimensionedScalar Cg4_;
            dimensionedScalar Cg5_;

            volScalarField KolmogorovFlag_;
    
      // Fields
            volSymmTensorField R_;
            volScalarField k_;
            volScalarField epsilon_;

       // Forcing fields
            volSymmTensorField QR_;
            volScalarField     Qepsilon_;
            Switch          directImpose_;
            labelList       lesCellLables_;
            scalarField     lesCellValuesK_;
            scalarField     lesCellValuesEpsilon_;
            symmTensorField lesCellValuesR_;
            label           imposeTurbEvery_;

        // Private member functions
        tmp<volTensorField> symm2full(volSymmTensorField& symm ) const;
        void updateKolmogorovFlag();
        tmp<volScalarField> T() const;
        tmp<volScalarField> L() const;

public:

    // Constructors

        //- Construct from components
        muGenElliptic
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const volSymmTensorField & RAvg,
            const volScalarField & epsilonAvg
        );


    //- Destructor
    virtual ~muGenElliptic()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return CmuKE_*sqr(k_)/(epsilon_ + epsilonSmall_);;
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
        virtual tmp<volSymmTensorField> R() const
        {
            return R_;
        }

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- update Kolmogorov flag field only
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
