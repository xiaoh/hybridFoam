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
#include "muRASModel.H"

// Added from channelFoam
#include "IFstream.H"
#include "OFstream.H"
#include "Random.H"
#include "scalar.H"
#include "chPressureGrad.H"
#include "autoPtr.H"
#include "ransForcing.H"
#include "Pstream.H"
#include "commChannel.H"
#include "splitRansFoamVersion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "printVersions.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createGradP.H"  // for channelFlow with constant mass flux only!

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    /* if(Pstream::master())
    {
        system("communicatorLR RANS > logRansComm &");
        } */
    Info<< "Starting time loop\n" << endl;
    
    // versionLRFoam defined in: "createFieldsLR.H"
    Info<< "Running splitRansFoam with " << versionLRFoam << "\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
      
        #include "readPISOControls.H"
        #include "CourantNo.H"
        
        // Pressure-velocity, combined PISO/SIMPLE corrector
        for (int oCorr=0; oCorr<nOuterCorr; oCorr++)
        {
            if (nOuterCorr != 1) { p.storePrevIter(); }
            // Momentum predictor
                #include "UEqn.H"  
            // --- PISO loop 
                #include "pEqn.H"
            ransModel->correct();
            if(RPtr.valid()) { RPtr() = ransModel->R()(); }
            gradP.adjust(rUA);
            
            // Finally correct relaxation forces
            relaxQ.correct();
        } // End of outer loop
        
        runTime.write();
        
                    
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
