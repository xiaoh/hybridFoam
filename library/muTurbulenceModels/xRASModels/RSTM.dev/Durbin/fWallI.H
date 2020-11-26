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

Global
    fWallI

Description
    Calculate wall f boundary conditions.

\*---------------------------------------------------------------------------*/

{
    labelList cellBoundaryFaceCount(f_.size(), 0);

    const fvPatchList& patches = mesh_.boundary();

    //- Initialise the near-wall f fields to zero
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isType<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                f_[faceCelli] = symmTensor::zero;
            }
        }
    }

    //- Accumulate the wall face contributions to f
    //  Increment cellBoundaryFaceCount for each face for averaging
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

         if (isType<wallFvPatch>(curPatch))
        {
            
            if ((f_.boundaryField()[patchi]).type() != "fixedInternalValue")
            {
                FatalErrorIn("boundary-condition evaluation")
                    << f_.boundaryField()[patchi].type()
                    << " is the wrong f patchField type for boundary-condition on patch "
                    << curPatch.name() << nl
                    << "    should be fixedInternalValue"
                    << exit(FatalError);
            }
            
            const scalarField& nuw = nu().boundaryField()[patchi];

            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];
                
                // For corner cells (with two boundary or more faces),
                // f in the near-wall cell is calculated
                // as an average
                
                cellBoundaryFaceCount[faceCelli]++;

                symmTensor fbc = -20.0 * sqr(nuw[facei]) * R_[faceCelli] / epsilon_[faceCelli]
                    /sqr(sqr(RASModel::y_[patchi][facei]));

    		    f_[faceCelli].yy() += fbc.yy();
    		    f_[faceCelli].xy() += fbc.xy();
    		    f_[faceCelli].yz() += fbc.yz();

                f_[faceCelli].xx() += -0.5*fbc.yy();
                f_[faceCelli].zz() += -0.5*fbc.yy();

                f_[faceCelli].xz() *= 0.0;
            }
        }
    }


    // Perform the averaging
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isType<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];
                f_[faceCelli] /= cellBoundaryFaceCount[faceCelli];
            }
        }
    }
}


// ************************************************************************* //
