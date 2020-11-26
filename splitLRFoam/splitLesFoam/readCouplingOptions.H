Info << "Reading dictionary relaxParameters\n" << endl;
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

dictionary couplingDict(relaxParameters.subDictPtr("couplingOptions"));
