    if (runTime.outputTime())
    {
        OFstream gradPLFile
        (
            runTime.path()/runTime.timeName()/"uniform"/"gradPL.raw"
        );

        OFstream gradPRFile
        (
            runTime.path()/runTime.timeName()/"uniform"/"gradPR.raw"
        );

        if(gradPLFile.good())
        {
            gradPLFile << gradPL << endl;
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Cannot open file "
                << runTime.path()/runTime.timeName()/"uniform"/"gradPL.raw"
                << exit(FatalError);
        };

        if(gradPRFile.good())
        {
            gradPRFile << gradPR << endl;
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Cannot open file "
                << runTime.path()/runTime.timeName()/"uniform"/"gradPR.raw"
                << exit(FatalError);
        };

    };
