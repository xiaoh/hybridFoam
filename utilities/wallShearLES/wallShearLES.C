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
    wallShearLES

Description
    Calculates and reports wall shear stress on specified wall, for
    the specified times.  Currently only intended for use in *spanwise
    averaged LES* cases (2D).
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "LESModel.H"
#include "nearWallDist.H"
#include "wallDist.H"
#include "OFstream.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"

#   include "createFields.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        const fvPatchList& patches = mesh.boundary();
        
        const word refPatchName(wallShearDict.lookup("refWall"));
        dimensionedScalar nu(transportProperties.lookup("nu"));

        label refPatchI = -1;


        OFstream wallShearFile
          (
           runTime.path()/runTime.timeName()/"wallShear.xy"
           );

        if( wallShearFile.good())
          {
            wallShearFile << "# x,  y, tau (tang. component of shear / rho, unit: m^2/s^2)" << endl;
          }
        else
          FatalErrorIn(args.executable())
            << "Cannot open file "
            << runTime.path()/runTime.timeName()/"wallShear.xy"
            << exit(FatalError);
        
        // Find the reference patch and ref. tau
        forAll(patches, patchi)
        {
            const fvPatch& currPatch = patches[patchi];

            if (isA<wallFvPatch>(currPatch) &&
                currPatch.name() == refPatchName)
            {
                Info<< "Processing patch " << currPatch.name()
                    << endl;


                Foam::Field<Foam::Vector<double> > forwardN =
                  (mesh.Sf().boundaryField()[patchi]/mesh.magSf().boundaryField()[patchi])
                  ^ normalZ;

                Tau.boundaryField()[patchi] =
                  nu.value()*(UMean.boundaryField()[patchi].snGrad())  & forwardN;
                Info << UMean.boundaryField()[patchi] << endl;
                const fvPatchScalarField & tauWall = 
                  Tau.boundaryField()[patchi];

                const fvsPatchField<Foam::Vector<double> >  & Cwall = 
                  mesh.Cf().boundaryField()[patchi];

                refPatchI = patchi;
                Info << "Writing to file wallShear.xy" << nl << endl;

                forAll(Cwall, faceI)
                  {
                    wallShearFile << Cwall[faceI][0] << ", "
                    << Cwall[faceI][1] << ", "
                    << tauWall[faceI] << endl;
                  }

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
              << "Patch named " << refPatchName 
              <<  " not found, " << " or its type is not WALL" << nl 
              << "Please specifiy in constant/wallShearDict "
              << "or check your boundary file" << abort(FatalError);
          }
    }
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
