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

Application
    myPlusLES (based on yPlusLES)

Description
    Calculates and reports yPlus for all wall patches, for the specified times.
    Also report utau for the normalization of $u$ and $y$.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "wallFvPatch.H"
#include "nearWallDist.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        fvMesh::readUpdateState state = mesh.readUpdate();

        wallDist y(mesh, true);
        // Wall distance
        if (timeI == 0 || state != fvMesh::UNCHANGED)
        {
            Info<< "Calculating wall distance\n" << endl;
            Info<< "Writing wall distance to field "
                << y.name() << nl << endl;
            y.write();
        }


        volScalarField yPlus
        (
            IOobject
            (
                "yPlus",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("yPlus", dimless, 0.0)
        );

        volScalarField yStar
        (
            IOobject
            (
                "yStar",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("yStar", dimless, 0.0)
        );

        volVectorField uStar
        (
            IOobject
            (
                "uStar",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("uStar", dimless, vector::zero)
        );

        volVectorField UMeanStar
        (
            IOobject
            (
                "UMeanStar",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("UMeanStar", dimless, vector::zero)
        );


        // Utau for normalization
        volScalarField uTau
        (
            IOobject
            (
                "uTau",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("uTau", dimVelocity, 0.0)
        );


        Info<< "Reading field U\n" << endl;
        volVectorField U
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );


        Info<< "Reading field UMean\n" << endl;

        volVectorField UMean
          (
           IOobject
           (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
            ),
           mesh
           );
        

#       include "createPhi.H"

        singlePhaseTransportModel laminarTransport(U, phi);

//         autoPtr<incompressible::LESModel> sgsModel
//         (
//             incompressible::LESModel::New(U, phi, laminarTransport)
//         );

        volScalarField::GeometricBoundaryField d = nearWallDist(mesh).y();
        volScalarField nuEff = laminarTransport.nu();

        const fvPatchList& patches = mesh.boundary();
        
        // U_tau for all the walls averaged. For computing the
        // velocity and length scale.  NOTE: not ideal for now, should
        // look for U_tau corresponding to the nearest wall instead.
        dimensionedScalar utauAvg("utauAvg", dimVelocity, 0);

        scalar nPatch = 0;
                
        forAll(patches, patchi)
        {
            const fvPatch& currPatch = patches[patchi];

            if (typeid(currPatch) == typeid(wallFvPatch))
            {
                yPlus.boundaryField()[patchi] =
                    d[patchi]
                   *sqrt
                    (
                        nuEff.boundaryField()[patchi]
                       *mag(U.boundaryField()[patchi].snGrad())
                    )
                   /nuEff.boundaryField()[patchi];
                const scalarField& Yp = yPlus.boundaryField()[patchi];

                uTau.boundaryField()[patchi] =
                  sqrt
                  (
                   nuEff.boundaryField()[patchi]
                   *mag(UMean.boundaryField()[patchi].snGrad())
                   );

                const fvPatchScalarField & utauWall = 
                  uTau.boundaryField()[patchi];
                
                dimensionedScalar utauTmp("utauTmp", dimVelocity, average(utauWall));
                
                utauAvg += utauTmp;
                nPatch ++;

                Info<< "Patch " << patchi
                    << " named " << currPatch.name()
                    << " y+ : min: " << min(Yp) << " max: " << max(Yp)
                    << " average: " << average(Yp) << nl << endl;
            }
        }

        utauAvg /= nPatch; 
        
        Info << "Avg. shear velocity is: "
             << utauAvg << "; # of walls: " << nPatch << "." << endl;
        
        yStar = y.y() * utauAvg / nuEff;
        uStar = U / utauAvg;
        UMeanStar = UMean / utauAvg;
        
        Info << "Writing yStar, uStar, and UMeanStar to corresonding fields." << endl;
        yStar.write();
        uStar.write();
        UMeanStar.write();

        Info<< "Writing yPlus to field "
            << yPlus.name() << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
