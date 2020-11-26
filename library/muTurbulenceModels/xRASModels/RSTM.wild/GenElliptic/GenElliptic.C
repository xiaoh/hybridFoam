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

\*---------------------------------------------------------------------------*/

#include "GenElliptic.H"
#include "wallFvPatch.H"
#include "gaussLaplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
  makeFvLaplacianTypeScheme(gaussLaplacianScheme, symmTensor, symmTensor)
}
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

GenElliptic::GenElliptic
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),
    mesh_(U.mesh()),
    CmuKE_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CmuKE",
            coeffDict_,
            0.09
        )
    ),
    Clrr1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clrr1",
            coeffDict_,
            1.22
        )
    ),
    Clrr2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clrr2",
            coeffDict_,
            0.6
        )
    ),



    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            1.9
        )
    ),
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.23
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.65
        )
    ),
    sigmaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK",
            coeffDict_,
            1.2
        )
    ),

    CL_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CL",
            coeffDict_,
            0.28
        )
    ),
    CEta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CEta",
            coeffDict_,
            80.0
        )
    ),
    yStarLim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "yStarLim",
            coeffDict_,
            5.0
        )
    ),

    SSG_(coeffDict_.lookupOrAddDefault<Switch>("SSG", false)),
    Cg1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cg1",
            coeffDict_,
            3.4
        )
    ),

    Cg1s_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cg1s",
            coeffDict_,
            1.8
        )
    ),

    Cg2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cg2",
            coeffDict_,
            4.2
        )
    ),

    Cg3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cg3",
            coeffDict_,
            0.8
        )
    ),

    Cg3s_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cg3s",
            coeffDict_,
            1.3
        )
    ),

    Cg4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cg4",
            coeffDict_,
            1.25
        )
    ),

    Cg5_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cg5",
            coeffDict_,
            0.4
        )
    ),

    yw_(mesh_),
    R_
    (
        IOobject
        (
            "R",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{
    //    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Cmpt, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Tensor<Cmpt>, PatchField, GeoMesh> >
GenElliptic::symm2full
(
    GeometricField<SymmTensor<Cmpt>, PatchField, GeoMesh>& symm
)
{
    typedef GeometricField<Tensor<Cmpt>, PatchField, GeoMesh> FieldType;
    tmp<FieldType> tFull
    (
        new FieldType
        (
           IOobject
           (
               symm.name() + "FullTensor",
               symm.mesh().time().timeName(),
               symm.mesh(),
               IOobject::NO_READ,
               IOobject::NO_WRITE
           ),
           symm.mesh(),
           dimensionedTensor("0", symm.dimensions(), Tensor<Cmpt>::zero),
           symm.boundaryField().types()
        )
    );

    FieldType& full = tFull();

    // manipulate components
    full.replace(Tensor<Cmpt>::XX, symm.component(SymmTensor<Cmpt>::XX));
    full.replace(Tensor<Cmpt>::YY, symm.component(SymmTensor<Cmpt>::YY));
    full.replace(Tensor<Cmpt>::ZZ, symm.component(SymmTensor<Cmpt>::ZZ));
    full.replace(Tensor<Cmpt>::XY, symm.component(SymmTensor<Cmpt>::XY));
    full.replace(Tensor<Cmpt>::YZ, symm.component(SymmTensor<Cmpt>::YZ));
    full.replace(Tensor<Cmpt>::XZ, symm.component(SymmTensor<Cmpt>::XZ));
    full.replace(Tensor<Cmpt>::ZX, symm.component(SymmTensor<Cmpt>::XZ));
    full.replace(Tensor<Cmpt>::YX, symm.component(SymmTensor<Cmpt>::XY));
    full.replace(Tensor<Cmpt>::ZY, symm.component(SymmTensor<Cmpt>::YZ));

    return tFull;
}


tmp<volScalarField> GenElliptic::T() const
{
    volScalarField yStar_=pow(CmuKE_,0.25)*sqrt(k_)*yw_/nu();
    return max
    (
        k_/(epsilon_ + epsilonSmall_),
        pos(yStarLim_ - yStar_)*6.0*sqrt(nu()/(epsilon_ + epsilonSmall_))
    );
}


tmp<volScalarField> GenElliptic::L() const
{
    volScalarField yStar_=pow(CmuKE_,0.25)*sqrt(k_)*yw_/nu();
    return CL_*max
    (
        pow(k_,1.5)/(epsilon_ + epsilonSmall_),
        pos(yStarLim_ - yStar_)*CEta_
      * pow(pow(nu(),3.0)/(epsilon_ + epsilonSmall_), 0.25)
    );
}


tmp<volSymmTensorField> GenElliptic::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            R_ - nu()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> GenElliptic::divDevReff(volVectorField& U) const
{
    return fvc::div(R_) - fvm::laplacian(nu(), U);
}


bool GenElliptic::read()
{
    if (RASModel::read())
    {
        CmuKE_.readIfPresent(coeffDict());
        Clrr1_.readIfPresent(coeffDict());
        Clrr2_.readIfPresent(coeffDict());
        Cg1_.readIfPresent(coeffDict());
        Cg1s_.readIfPresent(coeffDict());
        Cg2_.readIfPresent(coeffDict());
        Cg3_.readIfPresent(coeffDict());
        Cg3s_.readIfPresent(coeffDict());
        Cg4_.readIfPresent(coeffDict());
        Cg5_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        Cmu_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        sigmaK_.readIfPresent(coeffDict());

        CEta_.readIfPresent(coeffDict());
        yStarLim_.readIfPresent(coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
