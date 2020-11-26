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
            
            vectorField wallNormalDir =
                (mesh_.Sf().boundaryField()[patchi]/mesh_.magSf().boundaryField()[patchi]);
            
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];
                
                cellBoundaryFaceCount[faceCelli]++;

                // wall normal direction is the y of local coordinate
                vector ydir = wallNormalDir[facei];
                vector xdir, zdir;

                if(wallsAlignedWithZ_)
                {
                    // assume spanwise direction of local coord is aligned with z axis of global coord
                    zdir = vector(0.0, 0.0, 1.0);
                    // use right-hand rule to obtain x basis
                    xdir = ydir ^ zdir;
                }
                else
                {
                    // find direction of velocity of the near-wall cell
                    xdir = U_[faceCelli];
                    if(mag(xdir) < SMALL)
                    {
                        Info << "muDurbin::fWallI::Warning: " << nl
                             << " I cannott find x direction for wall-aligned frame. Use global x instead." << endl;
                        xdir = vector(1, 0, 0);
                    }
                    xdir -= (xdir & ydir) * ydir ; // subtract wall-normal component
                    xdir /= mag(xdir) + VSMALL; // normalization
                    zdir = xdir ^ ydir;
                }

                // form coordinate transformation matrix
                // Theoretical notes:
                // transformation tensor: Qij = e_i . e'_j;  where:
                // e_i and e'_j are basis vectors in the two coordinate systems
                // e_i is I in our case; therefore, Qt should be as follows:
                tensor Qt
                    (
                        xdir.x(), ydir.x(), zdir.x(),
                        xdir.y(), ydir.y(), zdir.y(),
                        xdir.z(), ydir.z(), zdir.z()
                    );
                    
                // transform R to local coordinate 
                // NOTE: product of two "matrices" corresponds to inner product of two tensors!
                // Same below in the transformation of "f" back to global coordinate
                symmTensor Rp = symm((Qt.T() & R_[faceCelli]) & Qt);
                
                // f--R relation valid for a few components
                symmTensor fpbc = -20.0 * sqr(nuw[facei]) * Rp / epsilon_[faceCelli]
                    /sqr(sqr(muRASModel::y_[patchi][facei]));
                
                // Hanjalic BC for f
                symmTensor fpH
                    (
                        -0.5*fpbc.yy(),    fpbc.xy(),       0.0,
                                           fpbc.yy(),       fpbc.yz(),
                                                           -0.5*fpbc.yy()
                    );

                // Durbin original BC for f
                symmTensor fpD
                    (
                        0.0,    fpbc.xy(),        0.0,
                                fpbc.yy(),   fpbc.yz(),
                                            -fpbc.yy()
                    );

                // f in local coordinate (f' or fprime)
                symmTensor fp = (1.0 - durFlag_) * fpH + durFlag_ * fpD;
                
                // transform fp back to global coordinate; add to f_ field;
                f_[faceCelli] += symm( (Qt & fp) & Qt.T());;
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
