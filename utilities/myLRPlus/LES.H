
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
                "yPlusL",
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
                "yStarL",
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
                "uStarL",
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
                "ULMeanStar",
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
                "uTauL",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("uTauL", dimVelocity, 0.0)
        );


        Info<< "Reading field UL\n" << endl;
        volVectorField U
        (
            IOobject
            (
                "UL",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );


        Info<< "Reading field ULMean\n" << endl;

        volVectorField UMean
          (
           IOobject
           (
            "ULMean",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
            ),
           mesh
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
                   /sgsModel->nu().boundaryField()[patchi];
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
        
        yStar = y.y() * utauAvg / sgsModel->nu();
        uStar = U / utauAvg;
        UMeanStar = UMean / utauAvg;
        
        Info << "Writing yStar, uStar, and UMeanStar to corresonding fields." << endl;
        yStar.write();
        uStar.write();
        UMeanStar.write();
    }
