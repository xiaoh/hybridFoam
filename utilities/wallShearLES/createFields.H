        // Utau for normalization
        volScalarField Tau
        (
            IOobject
            (
                "Tau",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Tau", dimVelocity*dimVelocity, 0.0)
        );

        Info<< "Reading field UMean\n" << endl;
        IOdictionary wallShearDict
          (
           IOobject
           (
            "wallShearDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
            )
           );

fileName constantDir = "./";
if (! isDir(runTime.path() + "/RANS" ))
{
    constantDir = "../../constant";
}

Info << "constantDit = " << constantDir << endl;

        IOdictionary transportProperties
          (
           IOobject
           (
            "transportProperties",
            runTime.constant(),
            constantDir,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
            )
           );
        
        
const word UName(wallShearDict.lookupOrAddDefault<word>("UName", "UMean"));
vector normalZ(wallShearDict.lookupOrAddDefault<vector>("normalZ", vector(0, 0, 1)));

        volVectorField UMean
          (
           IOobject
           (
            UName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
            ),
           mesh
           );


