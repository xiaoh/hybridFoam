{
    label partnerFreq = 1;
    
    if(host == "RANS")
    {
        partnerFreq = couplingDict.lookupOrDefault("mapR2LEvery", 1);
    }
    else if(host == "LES")
    {
        partnerFreq = couplingDict.lookupOrDefault("mapL2REvery", 1);
    }

    IFstream str(runTime.path()/source/"system/controlDict"); 
    dictionary partnerCtrlDict(str);
    
    scalar pdt( partnerCtrlDict.lookupOrDefault<scalar>("deltaT", 0) );
    
    word myStride = runTime.timeName( fetchFreq * runTime.deltaT().value() );
    word partnerStride = runTime.timeName(partnerFreq * pdt);
    if( partnerStride != myStride )
    {
        FatalErrorIn(args.executable())
            << " The interval of communication of LES/RANS should be consistent!"
            << " You provided " << nl
            << " my stride  (deltaT X commEvery) : " << runTime.deltaT().value() << " X " << fetchFreq << nl
            << " partner stride                  : " << pdt << " X " << partnerFreq << nl
            << exit(FatalError);        
    }
    else
    {
        Info << "Common stride confirmed: " << myStride << endl;
    }
}
