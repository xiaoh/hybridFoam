/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "homogeneousSpanIndex.H"
#include "boolList.H"
#include "syncTools.H"
#include "OFstream.H"
#include "meshTools.H"
#include "Time.H"
#include "SortableList.H"
#include "surfaceFields.H"
// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculate regions.
void Foam::homogeneousSpanIndex::calcLayeredRegions
(
    const fvMesh& mesh
)
{
    boolList blockedFace(mesh.nFaces(), false);

    // Determines face blocking
    surfaceVectorField normals = mesh.Sf()/mesh.magSf();

    label nfInternal = mesh.nInternalFaces();
    forAll(mesh.faces(), faceI)
      {
        if(faceI >= nfInternal) break;
        vector sn = normals[faceI];
        blockedFace[faceI] = ( mag(sn & spanDirection_) < 1e-5);
      }

    Info << "Finished checking blockFace set." << endl;
    
    if (false)
    {
        OFstream str(mesh.time().path()/"blockedFaces.obj");
        label vertI = 0;
        forAll(blockedFace, faceI)
        {
            if (blockedFace[faceI])
            {
                const face& f = mesh.faces()[faceI];
                forAll(f, fp)
                {
                    meshTools::writeOBJ(str, mesh.points()[f[fp]]);
                }
                str<< 'f';
                forAll(f, fp)
                {
                    str << ' ' << vertI+fp+1;
                }
                str << nl;
                vertI += f.size();
            }
        }
    }


    // Do analysis for connected regions
    cellRegion_.reset(new regionSplit(mesh, blockedFace));

    Info<< "Detected " << cellRegion_().nRegions() << " layers." << nl << endl;

    // Sum number of entries per region
    regionCount_ = regionSum(scalarField(mesh.nCells(), 1.0));

    // Average cell centres to determine ordering.
    pointField regionCc
    (
        regionSum(mesh.cellCentres())
      / regionCount_
    );

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::homogeneousSpanIndex::homogeneousSpanIndex
(
    const fvMesh& mesh
)
:
spanDirection_(0, 0, 1)
{
    // Calculate regions.
    calcLayeredRegions(mesh);
}


// ************************************************************************* //
