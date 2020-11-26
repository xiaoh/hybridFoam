Info << "Reading dictionary relaxParameters\n" << endl;

if(host == "RANS")
{
    localDir = "../../LES/constant";
}
IOdictionary relaxParameters
(
    IOobject
    (
        "relaxParameters",
        runTime.constant(),
        localDir,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dictionary couplingDict(relaxParameters.subDictPtr("couplingOptions"));
