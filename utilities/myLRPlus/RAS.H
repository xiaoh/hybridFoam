
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
                "yPlusR",
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
                "yStarR",
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
                "uStarR",
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
                "URMeanStar",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("URMeanStar", dimless, vector::zero)
        );


        // Utau for normalization
        volScalarField uTau
        (
            IOobject
            (
                "uTauR",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("uTauR", dimVelocity, 0.0)
        );


        Info<< "Reading field UR\n" << endl;
        volVectorField U
        (
            IOobject
            (
                "UR",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );


        Info<< "Reading field URMean\n" << endl;

        volVectorField UMean
          (
           IOobject
           (
            "URMean",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
            ),
           mesh
           );
        
        #undef  createPhi_H
        #include "createPhi.H"

        singlePhaseTransportModel laminarTransport(U, phi);

        volSymmTensorField RAvg
          (
           IOobject
           (
            "RAvg",
            runTime.timeName(),
            "consistent",
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
            ),
           mesh,
           dimensionedSymmTensor
           ("RAvg", dimVelocity*dimVelocity, symmTensor::zero)
           );
        
        volScalarField epsilonAvg
          (
           IOobject
           (
            "epsilonAvg",
            runTime.timeName(),
            "consistent",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
            ),
           mesh,
           dimensionedScalar
           ("epsilonAvg", sqr(dimLength)/pow3(dimTime), scalar(0.0))
           );


        volScalarField lesFlagK
          (
           IOobject
           (
            "lesFlagK",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
            ),
           mesh,
           dimensionedScalar
           ("lesFlagK", dimless, scalar(1))
           );
        
        autoPtr<incompressible::muRASModel> rasModel
          (
           incompressible::muRASModel::New(U, phi, laminarTransport,
                                         RAvg, epsilonAvg)
           );

        volScalarField::GeometricBoundaryField d = nearWallDist(mesh).y();
        volScalarField nuEff = rasModel->nuEff();

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
                   /rasModel->nu().boundaryField()[patchi];
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
        
        yStar = y.y() * utauAvg / rasModel->nu();
        uStar = U / utauAvg;
        UMeanStar = UMean / utauAvg;
        
        Info << "Writing yStar, uStar, and UMeanStar to corresonding fields." << endl;
        yStar.write();
        uStar.write();
        UMeanStar.write();
    }
