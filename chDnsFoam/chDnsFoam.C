/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    chDnsFoam

Description
    Incompressible single-phase DNS solver for flow in a channel with constant
    mass flow rate.

Note
    The code of the solver is almost the same of the OpenFOAM standard solver 
    oodles, where the LES model was removed.

    Original author: Alberto Passalacqua
    Adopted for OF 1.6, Heng Xiao
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IFstream.H"
#include "emergencyExit.H"
#include "OFstream.H"
#include "Random.H"
// #include "Probe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "readTransportProperties.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"
#   include "createGradP.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    for(runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readPISOControls.H"
#       include "CourantNo.H"

	// Momentum predictor 
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
         ==
            flowDirection*gradP		
        );

        if (momentumPredictor)
        {
            solve(UEqn == -fvc::grad(p));
        }


        // --- PISO loop

        volScalarField rUA = 1.0/UEqn.A();

        for (int corr=0; corr<nCorr; corr++)
        {
            U = rUA*UEqn.H();
            phi = (fvc::interpolate(U) & mesh.Sf()) 
                + fvc::ddtPhiCorr(rUA, U, phi);

            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rUA, p) == fvc::div(phi)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (nonOrth == nNonOrthCorr)
                {
                    phi -= pEqn.flux();
                }
            }

#           include "continuityErrs.H"

            U -= rUA*fvc::grad(p);
            U.correctBoundaryConditions();
        }


        // Correct driving force for a constant mass flow rate

        // Extract the velocity in the flow direction
        dimensionedScalar magUbarStar = 
            (flowDirection & U)().weightedAverage(mesh.V());

        // Calculate the pressure gradient increment needed to 
        // adjust the average flow-rate to the correct value
        dimensionedScalar gragPplus = 
            (magUbar - magUbarStar)/rUA.weightedAverage(mesh.V());

        U += flowDirection*rUA*gragPplus;

        gradP += gragPplus;

        Info<< "Uncorrected Ubar = " << magUbarStar.value() << tab
            << "pressure gradient = " << gradP.value() << endl;

        runTime.write();

#       include "writeGradP.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        emergencyExitCase.execute();        
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
