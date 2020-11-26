        fvVectorMatrix UEqnL
          (
           fvm::ddt(UL)
           + fvm::div(phiL, UL)
           + sgsModel->divDevReff(UL)
           ==
           flowDirection*gradPL
           );
        
        fvVectorMatrix UEqnR
          (
           fvm::ddt(UR)
           + fvm::div(phiR, UR)
           + ransModel->divDevReff(UR)
           ==
           flowDirection*gradPR
           );

        UEqnL.relax();
        UEqnR.relax();
        
        if (momentumPredictor)
          {
            solve(UEqnL == -fvc::grad(pL));
            solve(UEqnR == -fvc::grad(pR));
          }



