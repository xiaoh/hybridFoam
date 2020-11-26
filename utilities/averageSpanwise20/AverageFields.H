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

#ifndef AverageFields_H
#define AverageFields_H

#include "GeometricField.H"
#include "IOobjectList.H"
#include "homogeneousSpanIndex.H"
#include "wordRe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
void AverageFields
(
    const fvMesh & meshSource,
    const fvMesh & meshTarget,
    const homogeneousSpanIndex & homogeneousSpanIndexing,
    const HashSet<word>& selectedFields
)
{

    IOobjectList objects(meshSource, meshSource.time().timeName());
    wordRe meanField(".*Mean.*", wordRe::REGEXP_NOCASE);

    word fieldClassName
    (
        GeometricField<Type, fvPatchField, volMesh>::typeName
    );

    IOobjectList fields = objects.lookupClass(fieldClassName);

    Info << nl << "Processing all " << fieldClassName << nl << endl;

    for
    (
        IOobjectList::iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        if
            (
                (  selectedFields.size() && selectedFields.found(fieldIter()->name()) )
                ||
                ( !selectedFields.size() && meanField.match( fieldIter()->name() )    )
                ||
                (  selectedFields.size() && selectedFields.found("all") )
            )
        {
            Info << "\t" << fieldIter()->name() << endl;
            IOobject fieldTargetIOobject
                (
                    fieldIter()->name() + "_2d",
                    meshTarget.time().timeName(),
                    meshTarget,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                );
            
            // Read field fieldSource
            GeometricField<Type, fvPatchField, volMesh> fieldSource
                (
                    *fieldIter(),
                    meshSource
                );
            
            // Read fieldTarget
            GeometricField<Type, fvPatchField, volMesh> fieldTarget
                (
                    fieldTargetIOobject,
                    meshTarget,
                    dimensioned<Type> 
                    (
                        "dummy",fieldSource.dimensions(), fieldSource[0]*0.0
                    ),
                    zeroGradientFvPatchField<Type>::typeName
                );
            
            // Average fields over the domain down to a plane
            Field<Type> fieldSrc2D = homogeneousSpanIndexing.collapse(fieldSource);
            fieldTarget.internalField() = fieldSrc2D;
            
            Info << "    Calculation boundary field assuming zeroGradient." << endl;
            fieldTarget.correctBoundaryConditions();

            // Write field
            fieldTarget.write();
        }
        else
        {
            // Info << "Skipped " << fieldIter()->name() << endl;            
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
