Info << "Reading dictionary relaxParameters\n" << endl;
IOdictionary relaxParameters
(
    IOobject
    (
        "relaxParameters",
        runTime.constant(),
        meshL,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dictionary couplingDict(relaxParameters.subDictPtr("couplingOptions"));

word solverOption
(
    couplingDict.lookupOrAddDefault<word>("solverOption", "hybrid")
);

Switch solveLES  (! (solverOption == "RANS") );
Switch solveRANS (! (solverOption == "LES") );
Switch singleMesh(couplingDict.lookupOrAddDefault<Switch>("singleMesh", false));

Info << "******* SOLVER OPTION ******" << nl
     << "LES  ON:         " << solveLES << nl
     << "RANS ON:         " << solveRANS << nl
     << "Same Mesh:       " << singleMesh << nl
     << "****************************"
     << nl << endl;
