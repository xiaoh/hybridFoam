    // Last update:
   string versionLRFoam("Time-stamp: <2012-09-28 23:41:08 xiaoh>");

    Info<< "LES: Reading field p\n" << endl;
    volScalarField pL
    (
        IOobject
        (
            "pL",
            runTime.timeName(),
            meshL,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        meshL
    );

    Info<< "RAS: Reading field p\n" << endl;
    volScalarField pR
    (
        IOobject
        (
            "pR",
            runTime.timeName(),
            meshR,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        meshR
    );


    Info<< "LES: Reading field U\n" << endl;
    volVectorField UL
    (
        IOobject
        (
            "UL",
            runTime.timeName(),
            meshL,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        meshL
    );

   
    Info<< "RAS: Reading/creating field U\n" << endl;
    volVectorField UR
    (
        IOobject
        (
            "UR",
            runTime.timeName(),
            meshR,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        meshR
    );

  volSymmTensorField B
    (
     IOobject
     (
      "B",
      runTime.timeName(),
      meshL,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     meshL,
     dimensionedSymmTensor
       ("B", dimVelocity*dimVelocity, symmTensor::zero)
     );
   
    #   include "createPhiLR.H"

label pRefCell = 0;
scalar pRefValue = 0.0;

// Read in and set values for pRefcell & pRefValue
// read from meshL, but used for both pL and pR (see pLEqn.H and pREqn.H)
// for consistency
setRefCell(pL, meshL.solutionDict().subDict("PISO"), pRefCell, pRefValue);

singlePhaseTransportModel laminarTransportL(UL, phiL);
singlePhaseTransportModel laminarTransportR(UR, phiR);


    #include "createRelaxFields.H"
    #include "createTurbulence.H"  // sgs/RANS model and relaxQ

// This part is for emergency exit only
emergencyExit emergencyExitCase(runTime);
