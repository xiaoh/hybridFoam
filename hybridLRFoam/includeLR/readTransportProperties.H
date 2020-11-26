    Info<< "Reading transportProperties\n" << endl;
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            meshL,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );


   dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );

