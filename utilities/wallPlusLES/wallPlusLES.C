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
#include "LESModel.H"
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
            "UMean",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
            ),
           mesh
           );
        

        IOdictionary wallPlusDict
          (
           IOobject
           (
            "wallPlusDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
            )
           );
        

#       include "createPhi.H"

        singlePhaseTransportModel laminarTransport(U, phi);

        autoPtr<incompressible::LESModel> sgsModel
        (
            incompressible::LESModel::New(U, phi, laminarTransport)
        );

        volScalarField::GeometricBoundaryField d = nearWallDist(mesh).y();
        volScalarField nuEff = sgsModel->nuEff();

        const fvPatchList& patches = mesh.boundary();
        
        // U_tau for all the walls averaged. For computing the
        // velocity and length scale.  NOTE: not ideal for now, should
        // look for U_tau corresponding to the nearest wall instead.
        dimensionedScalar utauRef("utauRef", dimVelocity, 0);

        volScalarField utauRefField
        (
            IOobject
            (
                "utauRefField",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("utauRefField", dimVelocity, 0.0)
        );

        const word refPatchName(wallPlusDict.lookup("refWall"));
        Switch writeStar(wallPlusDict.lookupOrAddDefault<Switch>("writeStarFields", false));

        label refPatchI = -1;

        // Find the reference patch and ref. Utau
        forAll(patches, patchi)
        {
            const fvPatch& currPatch = patches[patchi];

            if (typeid(currPatch) == typeid(wallFvPatch) &&
                currPatch.name() == refPatchName)
            {
                Info<< "Patch " << patchi
                    << " named " << currPatch.name()
                    << " is the reference wall." << nl << endl;

                uTau.boundaryField()[patchi] =
                  sqrt
                  (
                   nuEff.boundaryField()[patchi]
                   *mag(UMean.boundaryField()[patchi].snGrad())
                   );

                const fvPatchScalarField & utauWall = 
                  uTau.boundaryField()[patchi];

                utauRef.value() = average(utauWall);

                refPatchI = patchi;
                break;
            }
        }

        if(refPatchI == -1)
          {
            Info << "Possible patches in the meshes are: "
                 << "[NAME (TYPE)]" << nl << endl;
            forAll(patches, patchi) {
              const fvPatch& currPatch = patches[patchi];
              Info << "Patch [" << patchi << "]: "
                   << currPatch.name() << " " 
                   << "(" << currPatch.type() 
                   << ")"  << endl;

              if (currPatch.name() == refPatchName)
                Info << "*** Patch with specified name found, but its type is not WALL!" << nl
                     << "Check file polyMesh/boundary. ***" << endl;
            }
            
            Info << nl << "Patch name you specified in the dictionary is: " 
                 << refPatchName << endl;

            FatalErrorIn("wallPlusLES")
              << "Reference patch named " << refPatchName 
              <<  " not found, " << " or its type is not WALL" << nl 
              << "Please specifiy in constant/wallPlusDict "
              << "or check your boundary file" << abort(FatalError);
          }


        scalar nPatch = 0;
        // Processing all the walls
        Info << "Processing all the walls ..." << endl;
        Info << "(using Utau from the reference wall only)" << endl;
        forAll(patches, patchi)
        {
          forAll (utauRefField.boundaryField()[patchi], celli) {
            utauRefField.boundaryField()[patchi][celli] = utauRef.value();
          }
          
          const fvPatch& currPatch = patches[patchi];
          
          if (typeid(currPatch) == typeid(wallFvPatch))
            {
              yPlus.boundaryField()[patchi] =
                d[patchi] * utauRefField.boundaryField()[patchi]
                / sgsModel->nu().boundaryField()[patchi];
              const scalarField& Yp = yPlus.boundaryField()[patchi];
              
              uTau.boundaryField()[patchi] =
                sqrt
                (
                 nuEff.boundaryField()[patchi]
                 *mag(UMean.boundaryField()[patchi].snGrad())
                 );
              
              const scalarField& ut = uTau.boundaryField()[patchi];
              
              nPatch ++;
              
              Info<< "Wall " << nPatch
                  << " (" << currPatch.name() << ")";
              if(patchi == refPatchI) Info << " [REFERENCE] ";
              Info << nl;
              Info << " y+   : min: " << min(Yp) << " max: " << max(Yp)
                   << " average: " << average(Yp) << nl
                   << " uTau : min: " << min(ut) << " max: " << max(ut)
                   << " average: " << average(ut) << nl << endl;
            }
        }
        
        Info << "Ref shear velocity: "
             << utauRef.value()  << ";  Shear stress: "
             << sqr(utauRef.value())
             << ";  # of walls: " << nPatch << "." << endl;

        if (writeStar) {
          yStar = y.y() * utauRef / sgsModel->nu();
          uStar = U / utauRef;
          UMeanStar = UMean / utauRef;
          
          Info << "Writing yStar, uStar, and UMeanStar to corresonding fields." << endl;
          yStar.write();
          uStar.write();
          UMeanStar.write();
          
          Info<< "Writing yPlus to field " << yPlus.name() << nl << endl;
          yPlus.write();
        }
    }
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
