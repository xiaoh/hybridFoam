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

#include "wallDistAndBoundBox.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace staticModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(wallDistAndBoundBox, 0);
addToRunTimeSelectionTable(staticModel, wallDistAndBoundBox, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

wallDistAndBoundBox::wallDistAndBoundBox
(
    const fvMesh& meshL,
    const fvMesh& meshR
)
:
    staticModel(typeName, meshL, meshR),

    ransRange_
        (
            "ransRange", 
            dimLength,
            paramDict_.lookup("ransRange")
        ),
    // by default: set min/max of bbox to zeros (and they overlap!)
    bboxLES_(),
    bboxRAS_()
{
    wallDist y(meshR_, true);
    lesFlagR_.internalField() = pos(y.y() - ransRange_);

    bool bblFound = paramDict_.readIfPresent("bboxLES", bboxLES_);
    bool bbrFound = paramDict_.readIfPresent("bboxRAS", bboxRAS_);

    if(bblFound && bbrFound && bboxLES_.overlaps(bboxRAS_))
    {
        Info << "Warning: specified LES and RAS bound boxes overlaps. " 
             << " Disable the bbox resolution model. Downgrade to wallVicinity model."
             << endl;
    }
    else
    {
        const  volVectorField & cellCenters = meshR_.C();
        if(meshR_.bounds().overlaps(bboxLES_) && bblFound)
        {
            forAll(cellCenters, celli)
            {
                if(bboxLES_.containsInside(cellCenters[celli]))
                {
                    lesFlagR_[celli] = 1.0;
                }
            }
        }
        
        if(meshR_.bounds().overlaps(bboxRAS_) && bbrFound)
        {
            forAll(cellCenters, celli)
            {
                if(bboxRAS_.containsInside(cellCenters[celli]))
                {
                    lesFlagR_[celli] = 0.0;
                }
            }
        }
    }

    interpLesFlagR2L();
    diagnosisFlagRegion();
    printParams();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool wallDistAndBoundBox::read()
{
    if (staticModel::read())
    {
        ransRange_.readIfPresent(paramDict());

        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace staticModels
} // End namespace Foam

// ************************************************************************* //
