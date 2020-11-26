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
#include "wallDistData.H"
#include "wallPointYPlus.H"
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
            0.25
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
    implicitDiv_(coeffDict_.lookupOrAddDefault<Switch>("implicitDiv", false)),

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

    KolmogorovFlag_
    (
        IOobject
        (
            "KolmogorovFlag",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("kolflag", dimless, 1.0)
    ),

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
// construct a tensor field out of a symmTensorField (by mirroring)
// dot product between two symmTensor is not correctly defined in OF.
Foam::tmp<Foam::volTensorField> GenElliptic::symm2full( volSymmTensorField& symm ) const
{
    tmp<volTensorField> tFull
        (  new volTensorField
           (
               IOobject
               (
                   symm.name()+"FullTensor",
                   runTime_.timeName(),
                   mesh_,
                   IOobject::NO_READ,
                   IOobject::NO_WRITE
               ),
               symm.mesh(),
               dimensionedTensor("DurbinEllitpic::ft", symm.dimensions(), tensor::zero),
               symm.boundaryField().types()
           )
        );

    volTensorField& full = tFull();

    // manipulate components
    full.replace(tensor::XX, symm.component(symmTensor::XX));
    full.replace(tensor::YY, symm.component(symmTensor::YY));
    full.replace(tensor::ZZ, symm.component(symmTensor::ZZ));
    full.replace(tensor::XY, symm.component(symmTensor::XY));
    full.replace(tensor::YZ, symm.component(symmTensor::YZ));
    full.replace(tensor::XZ, symm.component(symmTensor::XZ));
    full.replace(tensor::ZX, symm.component(symmTensor::XZ));
    full.replace(tensor::YX, symm.component(symmTensor::XY));
    full.replace(tensor::ZY, symm.component(symmTensor::YZ));

    return tFull;
}


void GenElliptic::updateKolmogorovFlag()
{
    // wall unit as defined by nu/sqrt(tauw/rho)
    volScalarField ystar
    (
        IOobject
        (
            "ystar",
            mesh_.time().constant(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("ystar", dimLength, 1.0)
    );

    const fvPatchList& patches = mesh_.boundary();
    forAll(patches, patchi)
    {
        if (isA<wallFvPatch>(patches[patchi]))
        {
            const fvPatchVectorField& Uw = U_.boundaryField()[patchi];
            const scalarField& nuw = nu().boundaryField()[patchi];
            // Note: nuw is used instead of nueff
            // for wall-resolving mesh, nut should be zero at wall
            ystar.boundaryField()[patchi] =
                nuw/sqrt(nuw*mag(Uw.snGrad()) + VSMALL);
        }
    }

    wallPointYPlus::yPlusCutOff = 500;
    wallDistData<wallPointYPlus> y(mesh_, ystar);
    
    KolmogorovFlag_ = pos(yStarLim_ - y/ystar);

    // For debug purpose only:
    if(runTime_.outputTime())
    {
        volScalarField KolmogorovFlagCompare
            (
                "KolmogorovFlagCompare", 
                pos(yStarLim_ - pow(CmuKE_,0.25)*sqrt(k_)*y/nu() )
            );
        KolmogorovFlagCompare.write();
    }

}

tmp<volScalarField> GenElliptic::T() const
{
    return max
        (
            k_/(epsilon_ + epsilonSmall_),
            KolmogorovFlag_ * 6.0 * sqrt(nu()/(epsilon_ + epsilonSmall_))
        );
}

tmp<volScalarField> GenElliptic::L() const
{
    return
        CL_*max
        (
            pow(k_,1.5)/(epsilon_ + epsilonSmall_),
            KolmogorovFlag_ * CEta_ * pow(pow(nu(),3.0)/(epsilon_ + epsilonSmall_),0.25)
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
    if(implicitDiv_)
    {
        return
            (
                fvc::div(R_)
                + fvc::laplacian(nut(), U, "laplacian(nuEff,U)")
                - fvm::laplacian(nuEff(), U)
            );
    }
    else
    {
        return
            (
                fvc::div(R_)
                - fvm::laplacian(nu(), U)
            );
    }
}


void GenElliptic::correct()
{
    updateKolmogorovFlag();
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
