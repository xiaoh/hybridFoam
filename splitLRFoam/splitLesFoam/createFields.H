    // Last update:
   string versionFoam("Time-stamp: <2012-01-02 11:46:50 xiaoh>");

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

    Info<< "LES: Reading field p\n" << endl;
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

    Info<< "LES: Reading field U\n" << endl;
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
   
  volSymmTensorField B
    (
     IOobject
     (
      "B",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh,
     dimensionedSymmTensor
       ("B", dimVelocity*dimVelocity, symmTensor::zero)
     );
   
#   include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;

// Read in and set values for pRefcell & pRefValue
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

singlePhaseTransportModel laminarTransport(U, phi);

#   include "readCouplingOptions.H"
#   include "createRelaxFields.H"
#   include "createTurbulence.H"  // sgs/RANS model and relaxQ

// This part is for emergency exit only
emergencyExit emergencyExitCase(runTime);
