  volScalarField cellDist
    (
     IOobject
     (
      "cellDist",
      "0",
      meshL,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
      ),
     meshL
     );


  labelIOList cellDecompositionR
    (
     IOobject
     (
      "cellDecompositionR",
      meshR.facesInstance(),
      meshR,
      IOobject::NO_READ,
      IOobject::NO_WRITE,
      false
      ),
     List<label>(meshR.cells().size(), -1)
     );

  // Write as volScalarField for postprocessing.
  volScalarField cellDistR
    (
     IOobject
     (
      "cellDistR",
      "0",
      meshR,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
     meshR,
     dimensionedScalar("cellDistR", dimless, -1),
     zeroGradientFvPatchScalarField::typeName
     );

  volScalarField procBorderFlag
    (
     IOobject
     (
      "procBorderFlag",
      "0",
      meshR,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
     meshR,
     dimensionedScalar("procBorderFlag", dimless, 0),
     zeroGradientFvPatchScalarField::typeName
     );
