/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

#include "muOmegaWallFunctionFvPatchScalarField.H"
#include "muRASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void muOmegaWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn("muOmegaWallFunctionFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

muOmegaWallFunctionFvPatchScalarField::muOmegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(p, iF),
    UName_("UR"),
    kName_("kR"),
    GName_("muRASModel::G"),
    nuName_("nu"),
    nutName_("nutR"),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    beta1_(0.075)
{
    checkType();
}


muOmegaWallFunctionFvPatchScalarField::muOmegaWallFunctionFvPatchScalarField
(
    const muOmegaWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedInternalValueFvPatchField<scalar>(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    kName_(ptf.kName_),
    GName_(ptf.GName_),
    nuName_(ptf.nuName_),
    nutName_(ptf.nutName_),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    beta1_(ptf.beta1_)
{
    checkType();
}


muOmegaWallFunctionFvPatchScalarField::muOmegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedInternalValueFvPatchField<scalar>(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "UR")),
    kName_(dict.lookupOrDefault<word>("k", "kR")),
    GName_(dict.lookupOrDefault<word>("G", "muRASModel::G")),
    nuName_(dict.lookupOrDefault<word>("nu", "nu")),
    nutName_(dict.lookupOrDefault<word>("nut", "nutR")),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    beta1_(dict.lookupOrDefault<scalar>("beta1", 0.075))
{
    checkType();
}


muOmegaWallFunctionFvPatchScalarField::muOmegaWallFunctionFvPatchScalarField
(
    const muOmegaWallFunctionFvPatchScalarField& owfpsf
)
:
    fixedInternalValueFvPatchField<scalar>(owfpsf),
    UName_(owfpsf.UName_),
    kName_(owfpsf.kName_),
    GName_(owfpsf.GName_),
    nuName_(owfpsf.nuName_),
    nutName_(owfpsf.nutName_),
    Cmu_(owfpsf.Cmu_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    beta1_(owfpsf.beta1_)
{
    checkType();
}


muOmegaWallFunctionFvPatchScalarField::muOmegaWallFunctionFvPatchScalarField
(
    const muOmegaWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(owfpsf, iF),
    UName_(owfpsf.UName_),
    kName_(owfpsf.kName_),
    GName_(owfpsf.GName_),
    nuName_(owfpsf.nuName_),
    nutName_(owfpsf.nutName_),
    Cmu_(owfpsf.Cmu_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    beta1_(owfpsf.beta1_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void muOmegaWallFunctionFvPatchScalarField::updateCoeffs()
{
    const muRASModel& rasModel = db().lookupObject<muRASModel>("RASProperties");
    const scalar yPlusLam = rasModel.yPlusLam(kappa_, E_);
    const scalarField& y = rasModel.y()[patch().index()];

    const scalar Cmu25 = pow(Cmu_, 0.25);

    volScalarField& G = const_cast<volScalarField&>
        (db().lookupObject<volScalarField>(GName_));

    volScalarField& omega = const_cast<volScalarField&>
        (db().lookupObject<volScalarField>(dimensionedInternalField().name()));

    const scalarField& k = db().lookupObject<volScalarField>(kName_);

    const scalarField& nuw =
        patch().lookupPatchField<volScalarField, scalar>(nuName_);

    const scalarField& nutw =
        patch().lookupPatchField<volScalarField, scalar>(nutName_);

    const fvPatchVectorField& Uw =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const scalarField magGradUw = mag(Uw.snGrad());

    // Set omega and G
    forAll(nutw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        scalar yPlus = Cmu25*y[faceI]*sqrt(k[faceCellI])/nuw[faceI];

        scalar omegaVis = 6.0*nuw[faceI]/(beta1_*sqr(y[faceI]));

        scalar omegaLog = sqrt(k[faceCellI])/(Cmu25*kappa_*y[faceI]);

        omega[faceCellI] = sqrt(sqr(omegaVis) + sqr(omegaLog));

        if (yPlus > yPlusLam)
        {
            G[faceCellI] =
                (nutw[faceI] + nuw[faceI])
               *magGradUw[faceI]
               *Cmu25*sqrt(k[faceCellI])
               /(kappa_*y[faceI]);
        }
        else
        {
            G[faceCellI] = 0.0;
        }
    }

    // TODO: perform averaging for cells sharing more than one boundary face
}


void muOmegaWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedInternalValueFvPatchField<scalar>::write(os);
    writeEntryIfDifferent<word>(os, "U", "UR", UName_);
    writeEntryIfDifferent<word>(os, "k", "kR", kName_);
    writeEntryIfDifferent<word>(os, "G", "muRASModel::G", GName_);
    writeEntryIfDifferent<word>(os, "nu", "nu", nuName_);
    writeEntryIfDifferent<word>(os, "nut", "nutR", nutName_);
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    os.writeKeyword("beta1") << beta1_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    muOmegaWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
