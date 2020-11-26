//
// createMesh.H
// ~~~~~~~~~~~~

    Foam::Info
        << "Create LES mesh for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;

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

    #include "readCouplingOptions.H"


autoPtr<Foam::fvMesh> meshRPtr(NULL);

if(!singleMesh)
{
    Foam::Info
        << "Create RANS mesh for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;

    meshRPtr.reset(
        new Foam::fvMesh
        (
            Foam::IOobject
            (
                "RANS",
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
            )
        )
    );
    if(!meshRPtr.valid())
        
        FatalErrorIn
        (
            "meshR initialization"
        )   << " creating RANS mesh failed"
            << endl
            << exit(FatalError);
}
else
{
    Info << "RANS mesh not created as singleMesh flag is on (by default)" << endl;
}

const Foam::fvMesh& meshR = singleMesh? meshL : meshRPtr();
