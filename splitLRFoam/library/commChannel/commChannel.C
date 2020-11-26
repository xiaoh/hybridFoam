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

  Class
  commChannel

  \*---------------------------------------------------------------------------*/

#include "commChannel.H"
#include "fvCFD.H"
#include "meshToMesh.H"
#include "MapVolFields.H"

//added for the reconstruction
#include "IOobjectList.H"
#include "processorMeshes.H"
#include "fvFieldReconstructor.H"

//added for the decomposition
#include "readFields.H"
#include "fvFieldDecomposer.H"


namespace Foam {

  // * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
  // * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

  // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
  commChannel::commChannel
  (
      Time& runTime,
      fileName& mappingSource,
      HashTable<word> fieldsToMap,
      HashSet<word> fieldsToReconstruct,
      HashSet<word> fieldsToDecompose
   )
  :
      runTime_(runTime),
      rootDirTarget_(runTime.rootPath()),
      caseDirTarget_(runTime.caseName().components()[0]),
      casePath_(mappingSource),
      rootDirSource_(casePath_.path()),
      caseDirSource_(casePath_.name()),
      fieldsToMap_(fieldsToMap),
      fieldsToReconstruct_(fieldsToReconstruct),
      fieldsToDecompose_(fieldsToDecompose),
      runTimeSource_
      (
          Time::controlDictName,
          rootDirSource_,
          caseDirSource_
      ),
      runTimeTarget_
      (
          Time::controlDictName,
          rootDirTarget_,
          caseDirTarget_
      )
  {
      Info<< "Source: " << rootDirSource_ << " " << caseDirSource_ << nl
          << "Target: " << rootDirTarget_ << " " << caseDirTarget_ << endl;
  }
  
  // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


  // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void commChannel::mapConsistentMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget
)
{
    // Create the interpolation scheme
    meshToMesh meshToMeshInterp
        (
            meshSource,
            meshTarget
        );

    Info << "Mapping fields for time " << meshSource.time().timeName()
         << endl;
    
    {
        IOobjectList objects(fieldsToMap_.size());
        wordList fieldList(fieldsToMap_.toc());
        forAll(fieldList, fieldI)
        {
            word fieldName = fieldList[fieldI];
            IOobject* objectPtr = new IOobject
                (
                    fieldName,
                    meshSource.time().timeName(), 
                    meshSource, 
                    IOobject::MUST_READ,  
                    IOobject::NO_WRITE
                );
            
            if (objectPtr->headerOk())
            {
                objects.insert(fieldName, objectPtr);
            }
            else
            {
                FatalErrorIn("commChannel::()")
                    << "Can not locate " << fieldName
                    << exit(FatalError);
            }
        }

        // Map volFields
        // ~~~~~~~~~~~~~
        MapVolFields<scalar>(objects, meshToMeshInterp, fieldsToMap_);
        MapVolFields<symmTensor>(objects, meshToMeshInterp, fieldsToMap_);
        MapVolFields<vector>(objects, meshToMeshInterp, fieldsToMap_);
    }

}


void commChannel::reconstruct()
{
    // determine the processor count directly
    label nProcs = 0;

    while (isDir(rootDirTarget_/caseDirTarget_/(word("processor") + name(nProcs))))
    {
        ++nProcs;
    }

    if (!nProcs)
    {
        FatalErrorIn("commChannel::reconstruct()")
            << "No processor* directories found"
            << exit(FatalError);
    }

    // Create the processor databases
    PtrList<Time> databases(nProcs);

    forAll (databases, procI)
    {
        databases.set
            (
                procI,
                new Time
                (
                    Time::controlDictName,
                    rootDirTarget_,
                    caseDirTarget_/fileName(word("processor") + name(procI))
                )
            );
    }

    // use the times list from the master processor
    // and select a subset based on the command-line options
    instant currentTime(runTime_.time().value());
    word regionName(fvMesh::defaultRegion);

    Foam::fvMesh mesh
        (
            Foam::IOobject
            (
                regionName,
                runTime_.timeName(),
                runTime_,
                Foam::IOobject::MUST_READ
            )
        );
    
    // Set all times on processor meshes equal to reconstructed mesh
    forAll (databases, procI)
    {
        databases[procI].setTime(runTime_.timeName(), runTime_.timeIndex());
    }

    // Read all meshes and addressing to reconstructed mesh
    processorMeshes procMeshes(databases, regionName);

    // Set time for global database
    //    runTime_.setTime(currentTime);
    
    Info << "\n(Reconstructor) Time = " << runTime_.timeName() 
         << endl;
    
    // Set time for all databases
    //    forAll (databases, procI)
    //    {
    //    databases[procI].setTime(currentTime);
    //    }
    
    // Get list of objects from processor0 database
    IOobjectList objects(procMeshes.meshes()[0], databases[0].timeName());

    // If there are any FV fields, reconstruct them
    if
        (
            objects.lookupClass(volScalarField::typeName).size()
            || objects.lookupClass(volVectorField::typeName).size()
            || objects.lookupClass(volSymmTensorField::typeName).size()
        )
    {
        Info << "Reconstructing FV fields: " << endl;
        
        fvFieldReconstructor fvReconstructor
            (
                mesh,
                procMeshes.meshes(),
                procMeshes.faceProcAddressing(),
                procMeshes.cellProcAddressing(),
                procMeshes.boundaryProcAddressing()
            );
        
        fvReconstructor.reconstructFvVolumeFields<scalar>
            (
                objects,
                fieldsToReconstruct_
            );
        fvReconstructor.reconstructFvVolumeFields<vector>
            (
                objects,
                fieldsToReconstruct_
            );
        fvReconstructor.reconstructFvVolumeFields<symmTensor>
            (
                objects,
                fieldsToReconstruct_
            );
    }
}


void commChannel::decompose()
{

    word regionName = fvMesh::defaultRegion;
    word regionDir = word::null;

    Info<< "\n(Decomposer) Time = " << runTime_.timeName() << endl;

    // determine the existing processor count directly
     label nProcs = 0;
     while
     (
         isDir
         (
             runTime_.path()
            /(word("processor") + name(nProcs))
            /runTime_.constant()
            /regionDir
            /polyMesh::meshSubDir
         )
     )
     {
         ++nProcs;
     }

    Info<< "Create mesh" << endl;
    fvMesh mesh
    (
        IOobject
        (
            regionName,
            runTime_.timeName(),
            runTime_
        )
    );

    // Search for list of objects for this time
    IOobjectList objects(mesh, runTime_.timeName());

    // Construct the vol fields
    // ~~~~~~~~~~~~~~~~~~~~~~~~
    PtrList<volScalarField> volScalarFields;
    readFields(mesh, objects, volScalarFields, fieldsToDecompose_);

    PtrList<volVectorField> volVectorFields;
    readFields(mesh, objects, volVectorFields, fieldsToDecompose_);

    PtrList<volSymmTensorField> volSymmTensorFields;
    readFields(mesh, objects, volSymmTensorFields, fieldsToDecompose_);

    // split the fields over processors
    for (label procI = 0; procI < nProcs; procI++)
    {
        Info<< "Processor " << procI << ": field transfer" << endl;

        // open the database
        Time processorDb
        (
            Time::controlDictName,
            runTime_.rootPath(),
            runTime_.caseName()/fileName(word("processor") + name(procI))
        );

        processorDb.setTime(runTime_);

        // read the mesh
        fvMesh procMesh
        (
            IOobject
            (
                regionName,
                processorDb.timeName(),
                processorDb
            )
        );

        labelIOList cellProcAddressing
        (
            IOobject
            (
                "cellProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        labelIOList boundaryProcAddressing
        (
            IOobject
            (
                "boundaryProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // FV fields
        if
        (
            volScalarFields.size()
         || volVectorFields.size()
         || volSymmTensorFields.size()
        )
        {
            labelIOList faceProcAddressing
            (
                IOobject
                (
                    "faceProcAddressing",
                    procMesh.facesInstance(),
                    procMesh.meshSubDir,
                    procMesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            fvFieldDecomposer fieldDecomposer
            (
                mesh,
                procMesh,
                faceProcAddressing,
                cellProcAddressing,
                boundaryProcAddressing
            );

            fieldDecomposer.decomposeFields(volScalarFields);
            fieldDecomposer.decomposeFields(volVectorFields);
            fieldDecomposer.decomposeFields(volSymmTensorFields);
        }

    }

}


// Get field from the other party (assumed to be already constructed)
// map to my mesh and write out (leave decompositio for later)
void commChannel::fetch(scalar sourceValue)
{ 
    instantList sourceTimes = runTimeSource_.times();
    label sourceTimeIndex = runTimeSource_.timeIndex();
    
    sourceTimeIndex = Time::findClosestTimeIndex
        (
            sourceTimes,
            sourceValue
        );
    
    if (sourceTimes[sourceTimeIndex].name() != runTime_.timeName())
    {

        FatalErrorIn("commChannel::fetch()")
            << "Looking for " << runTime_.timeName()
            << ", found only " << sourceTimes[sourceTimeIndex].name()  << endl
            << " I will quit without action now" 
            << exit(FatalError);
        
        return;
    }

    runTimeSource_.setTime(sourceTimes[sourceTimeIndex], sourceTimeIndex);
    runTimeTarget_.setTime(runTimeSource_);
    
    Info<< "\n(MapField) Source time: " << runTimeSource_.value()
        << ";  Target time: " << runTimeTarget_.value()
        << endl;

    Info<< "Create meshes (target)...  " ;
    fvMesh meshTarget
        (
            IOobject
            (
                fvMesh::defaultRegion,
                runTimeTarget_.timeName(),
                runTimeTarget_
            )
        );

    Info << "Create meshes (source)..." << endl;
    // Should be made class member as well.
    fvMesh meshSource
        (
            IOobject
            (
                fvMesh::defaultRegion,
                runTimeSource_.timeName(),
                runTimeSource_
            )
        );
    
    Info<< "Source mesh size: " << meshSource.nCells() << tab
        << "Target mesh size: " << meshTarget.nCells() << endl;
    
    mapConsistentMesh(meshSource, meshTarget);
}


  // * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


  // * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


  // * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace Foam


// ************************************************************************* //
