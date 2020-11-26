Info << "Reading dictionary relaxParameters from LES\n" << endl;
IOdictionary relaxParameters
(
    IOobject
    (
        "relaxParameters",
        runTime.constant(),
        "../../LES/constant",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dictionary couplingDict(relaxParameters.subDictPtr("couplingOptions"));
