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
    hybridLRFoam

Description
    Transient solver for incompressible flow with hybrid LES/RANS capability.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "muLESModel.H"
#include "muRASModel.H"
#include "emergencyExit.H"

// for mesh-to-mesh interpolation
#include "meshToMesh.H"
#include "meshTalkChannel.H"
#include "processorFvPatch.H"

// Added from channelFoam
#include "IFstream.H"
#include "OFstream.H"
#include "Random.H"
#include "scalar.H"
#include "chPressureGrad.H"
#include "autoPtr.H"
#include "relaxForcing.H"
#include "hybridLRFoamVersion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "printVersions.H"
    #include "createTime.H"
    #include "createMeshLR.H"
    #include "readTransportProperties.H" // Read nu & also Ubar
    #include "createFieldsLR.H"
    #include "initContinuityErrsLR.H"
    #include "createGradP.H"  // for channelFlow with constant mass flux only!
    #include "createConvectionBlending.H" // blended central/upwind for convection term in LES equation

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
   
    Info<< "Starting time loop\n" << endl;
    
    // versionLRFoam defined in: "createFieldsLR.H"
    //     Info<< "Running hybridLRFoam with " << versionLRFoam << "\n" << endl;

    scalarList splitTime(3, 0.0);

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
      
            #include "readPISOControlsLR.H"
            #include "CourantNoLR.H"
        
        runTime.cpuTimeIncrement();
        // Pressure-velocity, combined PISO/SIMPLE corrector
        for (int oCorr=0; oCorr<nOuterCorr; oCorr++)
        {
            if(solveLES) 
            {
                if (nOuterCorr != 1) { pL.storePrevIter(); }
                // Momentum predictor
                    #include "ULEqn.H"  
                // --- PISO loop 
                    #include "pLEqn.H"
                sgsModel->correct(); 
                B = sgsModel->B()();
                // Correct driving force for a constant mass flow rate
                gradPL.adjust(rUAL);
                Info <<"->LES done<-"<< endl;
            }
            splitTime[0] += runTime.cpuTimeIncrement();

            if(solveRANS)
            {
                if (nOuterCorr != 1) { pR.storePrevIter(); }
                // Momentum predictor
                    #include "UREqn.H"  
                // --- PISO loop 
                    #include "pREqn.H"
                ransModel->correct();
                if(RPtr.valid()) { RPtr() = ransModel->R()(); }
                gradPR.adjust(rUAR);
                Info <<"->RANS done<-"<< endl;
            }
            splitTime[1] += runTime.cpuTimeIncrement();

            // Finally correct relaxation forces
            relaxQ.correct();
            splitTime[2] += runTime.cpuTimeIncrement();
        } // End of outer loop
        
        runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << "  LES/RANS/Forcing = (" << splitTime[0] << ", "
            << splitTime[1] << ", "
            << splitTime[2] << ") s"
            << nl << endl;
        
        emergencyExitCase.execute();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
