  Info  << "Creating fields related to relaxation forcing" << endl;
  volVectorField UAvg
    (
     IOobject
     (
      "UAvg",
      runTime.timeName(),
      meshL,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
      ),
     UL
     );

    volSymmTensorField RResolvedAvg
    (
     IOobject
     (
      "RResolvedAvg",
      runTime.timeName(),
      meshL,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
      ),
     meshL,
     dimensionedSymmTensor
     ("RResolvedAvg", dimVelocity*dimVelocity, symmTensor::zero)
     );

    volSymmTensorField RSgsAvg
    (
     IOobject
     (
      "RSgsAvg",
      runTime.timeName(),
      meshL,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
      ),
     meshL,
     dimensionedSymmTensor
     ("RSgsAvg", dimVelocity*dimVelocity, symmTensor::zero)
     );

    volSymmTensorField RAvg
    (
     IOobject
     (
      "RAvg",
      runTime.timeName(),
      meshL,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
      ),
     RResolvedAvg + RSgsAvg
     );

    volScalarField epsilonAvg
    (
     IOobject
     (
      "epsilonAvg",
      runTime.timeName(),
      meshL,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
      ),
     meshL,
     dimensionedScalar
     ("epsilonAvg", sqr(dimLength)/pow3(dimTime), scalar(0))
     );


// Coupling force for LES
    volVectorField QLes
    (
     IOobject
     (
      "QLes",
       runTime.timeName(),
      meshL,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
     meshL,
     dimensionedVector("QLes", dimAcceleration, vector::zero)
     );

// Coupling force for Rans
    volVectorField QRans
    (
     IOobject
     (
      "QRans",
      runTime.timeName(),
      meshR,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     meshR,
     dimensionedVector("QRans", dimAcceleration, vector::zero)
     );

   // NO read on the mirror fields, since they will be updated before use anyway.
    volScalarField    mrEpsilonAvg
    (
        IOobject
        (
            "mrEpsilonAvg",
            runTime.timeName(),
            meshR,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        meshR,
        dimensionedScalar
        ("mrEpsilonAvg", sqr(dimLength)/pow3(dimTime), scalar(0.0))
    );

   volSymmTensorField mrRAvg
    (
        IOobject
        (
            "mrRAvg",
            runTime.timeName(),
            meshR,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        meshR,
        dimensionedSymmTensor
        ("mrRAvg", sqr(dimVelocity), symmTensor::zero)
    );

if(singleMesh)
{
    mrRAvg.writeOpt() = IOobject::NO_WRITE;
    mrEpsilonAvg.writeOpt() = IOobject::NO_WRITE;
}
