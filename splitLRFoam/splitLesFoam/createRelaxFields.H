  Info  << "Creating fields related to relaxation forcing" << endl;
  volVectorField UAvg
    (
     IOobject
     (
      "UAvg",
      runTime.timeName(),
      mesh,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
      ),
     U
     );

    volSymmTensorField RResolvedAvg
    (
     IOobject
     (
      "RResolvedAvg",
      runTime.timeName(),
      mesh,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
      ),
     mesh,
     dimensionedSymmTensor
     ("RResolvedAvg", dimVelocity*dimVelocity, symmTensor::zero)
     );

    volSymmTensorField RSgsAvg
    (
     IOobject
     (
      "RSgsAvg",
      runTime.timeName(),
      mesh,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
      ),
     mesh,
     dimensionedSymmTensor
     ("RSgsAvg", dimVelocity*dimVelocity, symmTensor::zero)
     );

    volSymmTensorField RAvg
    (
     IOobject
     (
      "RAvg",
      runTime.timeName(),
      mesh,
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
      mesh,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
      ),
     mesh,
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
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("QLes", dimAcceleration, vector::zero)
);

