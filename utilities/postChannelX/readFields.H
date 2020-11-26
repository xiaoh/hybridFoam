    IOobject UMeanHeader
    (
        "UMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (!UMeanHeader.headerOk())
    {
        Info<< "    No UMean field" << endl;
        continue;
    }

    volVectorField UMean
    (
        UMeanHeader,
        mesh
    );

    volSymmTensorField UPrime2Mean
    (
        IOobject
        (
            "UPrime2Mean",
            runTime.timeName(),
            mesh,
            // IOobject::MUST_READ
            IOobject::READ_IF_PRESENT
        ),
        mesh,
        symmTensor::zero
    );

    volScalarField Rxx(UPrime2Mean.component(symmTensor::XX));
    volScalarField Ryy(UPrime2Mean.component(symmTensor::YY));
    volScalarField Rzz(UPrime2Mean.component(symmTensor::ZZ));
    volScalarField Rxy(UPrime2Mean.component(symmTensor::XY));

    volScalarField pPrime2Mean
    (
        IOobject
        (
            "pPrime2Mean",
            runTime.timeName(),
            mesh,
            //            IOobject::MUST_READ
            IOobject::READ_IF_PRESENT
        ),
        mesh,
        0.0
    );


    volSymmTensorField BMean // SGS stress
    (
        IOobject
        (
            "BMean",
            runTime.timeName(),
            mesh,
            // IOobject::MUST_READ
            IOobject::READ_IF_PRESENT
        ),
        mesh,
        symmTensor::zero
    );

    volScalarField Bxx(BMean.component(symmTensor::XX));
    volScalarField Byy(BMean.component(symmTensor::YY));
    volScalarField Bzz(BMean.component(symmTensor::ZZ));
    volScalarField Bxy(BMean.component(symmTensor::XY));
