//
// ~~~~~~~~~~~~

    volVectorField U2
    (
        IOobject
        (
            "U2",
            runTime.timeName(),
            mesh2D,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh2D,
        dimensionedVector
        ("u2",  dimensionSet(0, 1, -1, 0, 0), vector::zero),
        fixedValueFvPatchField<vector>::typeName
    );

    volScalarField pPrime2Mean2
    (
     IOobject
     (
      "pPrime2Mean2",
      runTime.timeName(),
      mesh2D,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh2D,
     dimensionedScalar
     ("pPrime2Mean2",  dimensionSet(0, 2, -2, 0, 0), 0.0),
     fixedValueFvPatchField<scalar>::typeName
     );


    volSymmTensorField RResolved2
    (
     IOobject
     (
      "RResolved2",
      runTime.timeName(),
      mesh2D,
      IOobject::READ_IF_PRESENT,
      IOobject::NO_WRITE
      ),
     mesh2D,
     dimensionedSymmTensor
     ("RResolved2", dimVelocity*dimVelocity, symmTensor::zero)
     );

    volSymmTensorField B2
    (
     IOobject
     (
      "B2",
      runTime.timeName(),
      mesh2D,
      IOobject::READ_IF_PRESENT,
      IOobject::NO_WRITE
      ),
     mesh2D,
     dimensionedSymmTensor
     ("B2", dimVelocity*dimVelocity, symmTensor::zero)
     );

