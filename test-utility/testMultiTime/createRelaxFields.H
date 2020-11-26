  Info << "\n"
       << "Creating fields related to relaxation forcing\n" << endl;
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
     UL
     );

 

    volSymmTensorField RResolvedAvg
    (
     IOobject
     (
      "RResolvedAvg",
      runTime.timeName(),
      "consistent",
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
      "consistent",
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
      "consistent",
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh,
     dimensionedSymmTensor
       ("RAvg", dimVelocity*dimVelocity, symmTensor::zero)
     );

    volScalarField epsilonAvg
    (
     IOobject
     (
      "epsilonAvg",
      runTime.timeName(),
      "consistent",
      mesh,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
      ),
     mesh,
     dimensionedScalar
     ("epsilonAvg", sqr(dimLength)/pow3(dimTime), scalar(0.0))
     );

   volScalarField turbLength
    (
     IOobject
     (
      "turbLength",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh,
     dimensionedScalar
     ("turbLengths", dimLength, scalar(0.0))
     );

    volScalarField lesFlagL
    (
     IOobject
     (
      "lesFlagL",
       runTime.timeName(),
      mesh,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
      ),
     mesh,
     dimensionedScalar
     ("lesFlagL", dimless, scalar(1))
     );

    volScalarField lesFlagK
    (
     IOobject
     (
      "lesFlagK",
      runTime.timeName(),
      mesh,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
      ),
     mesh,
     dimensionedScalar
     ("lesFlagK", dimless, scalar(1))
     );

     volScalarField lesFlagGeneral
    (
     IOobject
     (
      "lesFlagGeneral",
      runTime.timeName(),
      mesh,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
      ),
     mesh,
     dimensionedScalar
     ("lesFlagGeneral", dimless, scalar(0))
     );

    volScalarField resIndicatorL
    (
     IOobject
     (
      "resIndicatorL",
       runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
     mesh,
     dimensionedScalar
     ("resIndicatorL", dimless, scalar(1))
     );

    volScalarField resIndicatorK
    (
     IOobject
     (
      "resIndicatorK",
       runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
     mesh,
     dimensionedScalar
     ("resIndicatorK", dimless, scalar(0))
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
     dimensionedVector
     ("QLes", dimAcceleration, vector::zero)
     );

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
     dimensionedVector
       (
        "QRans", dimAcceleration, vector::zero
       )
     );

   Info<< "Reading dictionary relaxParameters\n" << endl;
    IOdictionary relaxParameters
    (
        IOobject
        (
            "relaxParameters",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

