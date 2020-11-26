//
// createMesh.H
// ~~~~~~~~~~~~

    Foam::Info
        << "Create LES mesh for time = "
        << runTime.timeName() << Foam::endl;

    Foam::fvMesh meshL
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    Info<< "LES mesh size: " << meshL.nCells() << nl << endl;

    Foam::Info
        << "Create RANS mesh for time = "
        << runTime.timeName() << Foam::endl;

    Foam::fvMesh meshR
    (
        Foam::IOobject
        (
            "RANS",
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

Info<< "RANS mesh size: " << meshR.nCells() << endl;

IOdictionary decompositionDict
(
 IOobject
 (
  "decomposeParDict",
  runTime.system(),
  runTime,
  IOobject::MUST_READ,
  IOobject::NO_WRITE
  )
 );

int nProcs(readInt(decompositionDict.lookup("numberOfSubdomains")));
