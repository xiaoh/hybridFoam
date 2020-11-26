        // Initialization of rUA; Updated at the end of the loop
        volScalarField rUAL = 1.0/UEqnL.A();
        volScalarField rUAR = 1.0/UEqnR.A();
        
        for (int corr=0; corr<nCorr; corr++)
          {
            // Estimate U, considering only the "transport" and "source",
            // NOT including pressure gradient.
            // Ref. (Jasak PhD thesis, pp. 145)
            UL = rUAL*UEqnL.H();
            UR = rUAR*UEqnR.H();

            surfaceScalarField phiQL = fvc::interpolate(rUAL) 
                                    * (fvc::interpolate(QLes) & mesh.Sf());
            surfaceScalarField phiQR = fvc::interpolate(rUAR) 
                                    * (fvc::interpolate(QRans) & mesh.Sf());

            forAll(pL.boundaryField(), patchi)
              {
                if (isA<zeroGradientFvPatchScalarField>(pL.boundaryField()[patchi]))
                  {
                    phiQL.boundaryField()[patchi] = 0.0;
                    phiQR.boundaryField()[patchi] = 0.0;
                  }
              }

            phiL = (fvc::interpolate(UL) & mesh.Sf())
              + fvc::ddtPhiCorr(rUAL, UL, phiL)
              + phiQL;
            phiR = (fvc::interpolate(UR) & mesh.Sf())
              + fvc::ddtPhiCorr(rUAR, UR, phiR)
              + phiQR;
            
            adjustPhi(phiL, UL, pL);
            adjustPhi(phiR, UR, pR);
            
            // Non-orthogonal pressure corrector loop
            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
              {
                // Pressure corrector
                
                fvScalarMatrix pEqnL
                  (
                   fvm::laplacian(rUAL, pL) == fvc::div(phiL)
                   );
               
                fvScalarMatrix pEqnR
                  (
                   fvm::laplacian(rUAR, pR) == fvc::div(phiR)
                   );
               
 
                pEqnL.setReference(pRefCell, pRefValue);
                pEqnR.setReference(pRefCell, pRefValue);
                
                if
                  (
                   corr == nCorr-1
                   && nonOrth == nNonOrthCorr
                   )
                  {
                    pEqnL.solve(mesh.solver("pFinal"));
                    pEqnR.solve(mesh.solver("pFinal"));
                  }
                else
                  {
                    pEqnL.solve();
                    pEqnR.solve();
                  }
                
                if (nonOrth == nNonOrthCorr)
                  {
                    phiL -= pEqnL.flux();
                    phiR -= pEqnR.flux();
                  }
              }
            
            #include "continuityErrsLR.H"
            
            UL += fvc::reconstruct(phiQL) - rUAL*fvc::grad(pL);
            UL.correctBoundaryConditions();
            rUAL = 1.0/UEqnL.A();

            UR += fvc::reconstruct(phiQR) - rUAR*fvc::grad(pR);
            UR.correctBoundaryConditions();
            rUAR = 1.0/UEqnR.A();
          }
            
