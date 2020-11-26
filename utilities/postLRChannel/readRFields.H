word URMeanName = channelDict.lookupOrAddDefault<word>
(
    "URMeanName", "URMean" 
);

word RMeanName = channelDict.lookupOrAddDefault<word>
(
    "RRMeanName", "RRMean" 
);


    IOobject UMeanHeader
    (
        URMeanName,
        runTime.timeName(),
        meshR,
        IOobject::MUST_READ
    );

    if (!UMeanHeader.headerOk())
    {
        Info<< "    No URMean field" << endl;
        continue;
    }

    volVectorField UMean
    (
        UMeanHeader,
        meshR
    );

    volSymmTensorField RMean
    (
        IOobject
        (
            RMeanName,
            runTime.timeName(),
            meshR,
            IOobject::READ_IF_PRESENT
        ),
        meshR,
        symmTensor::zero
    );

    volScalarField RRxx(RMean.component(symmTensor::XX));
    volScalarField RRyy(RMean.component(symmTensor::YY));
    volScalarField RRzz(RMean.component(symmTensor::ZZ));
    volScalarField RRxy(RMean.component(symmTensor::XY));

    volScalarField kRMean
    (
        IOobject
        (
            "kR",
            runTime.timeName(),
            meshR,
            IOobject::MUST_READ
        ),
        meshR
    );

