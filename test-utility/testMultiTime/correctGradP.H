        // Extract the velocity in the flow direction
        dimensionedScalar magUbarStarL =
            (flowDirection & UL)().weightedAverage(mesh.V());

        // Calculate the pressure gradient increment needed to
        // adjust the average flow-rate to the correct value;
        // magUbar: target; magUbarStar: actual in simulation.
        dimensionedScalar gragPplusL =
            (magUbar - magUbarStarL)/rUAL.weightedAverage(mesh.V());

        UL += flowDirection*rUAL*gragPplusL;

        gradPL += gragPplusL;

        //------------------------------------------------------------
        // Extract the velocity in the flow direction
        dimensionedScalar magUbarStarR =
            (flowDirection & UR)().weightedAverage(mesh.V());

        // Calculate the pressure gradient increment needed to
        // adjust the average flow-rate to the correct value;
        // magUbar: target; magUbarStar: actual in simulation.
        dimensionedScalar gragPplusR =
            (magUbar - magUbarStarR)/rUAR.weightedAverage(mesh.V());

        UR += flowDirection*rUAR*gragPplusR;

        gradPR += gragPplusR;

       //------------------------------------------------------------
        Info<< "LES Uncorrected Ubar = " << magUbarStarL.value() << tab
            << "pressure gradient = " << gradPL.value() << endl;

        Info<< "RAN Uncorrected Ubar = " << magUbarStarR.value() << tab
            << "pressure gradient = " << gradPR.value() << endl;
