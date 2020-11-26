  Info  << "Creating fields related to relaxation forcing" << endl;

// Coupling force for Rans
    volVectorField QRans
    (
     IOobject
     (
      "QRans",
      runTime.timeName(),
      mesh,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedVector("QRans", dimAcceleration, vector::zero)
     );

   // NO read on the mirror fields, since they will be updated before use anyway.
    volScalarField    mrEpsilonAvg
    (
        IOobject
        (
            "mrEpsilonAvg",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        ("mrEpsilonAvg", sqr(dimLength)/pow3(dimTime), scalar(0.0))
    );

   volSymmTensorField mrRAvg
    (
        IOobject
        (
            "mrRAvg",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor
        ("mrRAvg", sqr(dimVelocity), symmTensor::zero)
    );
