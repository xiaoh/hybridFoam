/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    communicatorLR

Description


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "commChannel.H"
#include "IFstream.H"
#include "communicatorLRVersion.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    //    argList::validArgs.clean();
    argList::validArgs.append("host solver");

    #   include "setRootCase.H"
    #   include "printVersions.H"
    #   include "createTime.H"

    HashTable<word> fieldsToReceive(3);
    HashSet<word> fieldsToOffer(3);
    HashSet<word> fieldsToDecompose(3);
    fileName source;
    fileName localDir;
    label fetchFreq = 1;

    word host(args.additionalArgs()[0]);

    #include "readCouplingOptions.H"

    if(host == "RANS")
    {
        source = "../LES";
        fieldsToReceive.insert("UAvg", "mrUAvg");
        fieldsToReceive.insert("epsilonAvg", "mrEpsilonAvg");
        fieldsToReceive.insert("RAvg", "mrRAvg");
        fieldsToOffer.insert("U");
        fieldsToOffer.insert("epsilonR");
        fieldsToOffer.insert("RR");
        fieldsToDecompose.insert("mrUAvg");
        fieldsToDecompose.insert("mrEpsilonAvg");
        fieldsToDecompose.insert("mrRAvg");
        fetchFreq = couplingDict.lookupOrDefault("mapL2REvery", 1);
        Info << "Host solver type: RANS" << endl;
    }
    else if(host == "LES")
    {
        source = "../RANS";
        fieldsToReceive.insert("U", "mlUR");
        fieldsToReceive.insert("epsilonR", "mlEpsilonR");
        fieldsToReceive.insert("RR", "mlRR");
        fieldsToOffer.insert("UAvg");
        fieldsToOffer.insert("epsilonAvg");
        fieldsToOffer.insert("RAvg");
        fieldsToDecompose.insert("mlUR");
        fieldsToDecompose.insert("mlEpsilonR");
        fieldsToDecompose.insert("mlRR");
        fetchFreq = couplingDict.lookupOrDefault("mapR2LEvery", 1);
        Info << "Host solver type: LES" << endl;
    }
    else
    {
        FatalErrorIn(args.executable())
            << "Must specify host solver (LES or RANS)"
            << exit(FatalError);
    }
    
    #include "checkStride.H"

    commChannel mappingLR(runTime, source, fieldsToReceive, fieldsToOffer, fieldsToDecompose);

    Switch cleanFootprint(couplingDict.lookupOrAddDefault<Switch>("cleanFootprint", true) );
    
    fileName sigFileRoot(runTime.path());

    word prevt = runTime.time().timeName();
    bool prevWasOutTime = false;
    bool hostPar = isDir(sigFileRoot/"processor0");
    while (runTime.loop())
    {
        if(runTime.timeIndex() % fetchFreq) // not an exchange time
        {
            continue;
        }

        Info << nl << "Communicator Time: " << runTime.time().value() << endl;
        fileName hostReady(sigFileRoot/"hostReady"+runTime.time().timeName());
        fileName commReady(sigFileRoot/"commReady"+runTime.time().timeName());
        fileName partnerReadyMe(sigFileRoot/"partnerReady"+runTime.time().timeName());
        fileName partnerReadyYou(sigFileRoot/source/"partnerReady"+runTime.time().timeName());

        word startt = runTime.time().timeName(runTime.time().startTime().value());
        fileName currTime = sigFileRoot/"processor0"/runTime.time().timeName();

        Info << "Waiting for host ready ";
        while (1)
        {
            if(isFile(hostReady))
            {
                Info << " Done! " << endl;;
                    
                if(isDir(currTime))
                {
                    Info << " Reconstruct my fields." << endl;
                    mappingLR.reconstruct();
                }
                else
                {
                    Info << "Skipped reconstruction for serial case." << endl;
                }

                system("touch " + partnerReadyMe);
                break;
            }
            else
            {
                Info << ".";
                Foam::sleep(2);
            }
        }
        
        Info << "Waiting for partner ready ";
        while(1)
        {
            if(isFile(partnerReadyYou))
            {
                Info << " Done. \n" << " Map fields from partner." << endl;

                mappingLR.fetch(runTime.time().value());
                
                if(isDir(currTime))
                { 
                    mappingLR.decompose();
                }
                else
                {
                    Info << "Skipped decomposing for serial case." << endl;
                }

                system("touch " + commReady);
                break;
            }
            else
            {
                Info << ".";
                Foam::sleep(1);
            }
        }
        
        // Remove previous time (communication files)
        if(cleanFootprint && prevt!= startt)
        {
            fileName partnerPrev(sigFileRoot/"partnerReady"+prevt);
            Foam::rm(partnerPrev);
            
            fileName commPrev(sigFileRoot/"commReady"+prevt);
            Foam::rm(commPrev);
            
            //if(hostPar || (! prevWasOutTime) )
            if(! prevWasOutTime )
            {
                // Remove dir if:
                // 1. the host is parallel; OR
                // 2. the host is serial, but the previous time is not an output time
                Foam::rmDir(sigFileRoot/prevt);
            }
        }

        prevt = runTime.time().timeName();
        prevWasOutTime = runTime.outputTime();

        Info<< "\nExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << "----------------------------------------"
            << endl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
