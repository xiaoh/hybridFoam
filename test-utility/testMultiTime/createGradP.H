    dimensionedScalar gradPL
    (
        "gradPL",
        dimensionSet(0, 1, -2, 0, 0),
        0.0
    );

    dimensionedScalar gradPR
    (
        "gradPR",
        dimensionSet(0, 1, -2, 0, 0),
        0.0
    );

    IFstream gradPLFile
    (
        runTime.path()/runTime.timeName()/"uniform"/"gradPL.raw"
    );

    if(gradPLFile.good())
    {
        gradPLFile >> gradPL;
        Info<< "LES: Reading average pressure gradient" <<endl
            << endl;
    }
    else
    {
        Info<< "LES: Initializing with 0 pressure gradient" <<endl
            << endl;
    };


    IFstream gradPRFile
    (
        runTime.path()/runTime.timeName()/"uniform"/"gradPR.raw"
    );

    if(gradPRFile.good())
    {
        gradPRFile >> gradPR;
        Info<< "RANS: Reading average pressure gradient" <<endl
            << endl;
    }
    else
    {
        Info<< "RANS: Initializing with 0 pressure gradient" <<endl
            << endl;
    };
