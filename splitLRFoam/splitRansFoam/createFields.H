    // Last update:
   string versionLRFoam("Time-stamp: <2012-01-08 19:08:21 xiaoh>");

    Info<< "Reading transportProperties\n" << endl;
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

   dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );


    Info<< "RAS: Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
   
    Info<< "RAS: Reading/creating field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #   include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;

// Read in and set values for pRefcell & pRefValue
// for consistency
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

singlePhaseTransportModel laminarTransport(U, phi);

    #   include "readCouplingOptions.H"
    #   include "createRelaxFields.H"
    #   include "createTurbulence.H"  // sgs/RANS model and relaxQ
