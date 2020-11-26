// p equation RANS part

for (int corr=0; corr<nCorr; corr++)
{
    UR = rUAR*UEqnR.H();
    
    forAll(pR.boundaryField(), patchi)
    {
        if (isA<zeroGradientFvPatchScalarField>(pR.boundaryField()[patchi]))
        {
            phiQR.boundaryField()[patchi] = 0.0;
        }
    }
    
    phiR = (fvc::interpolate(UR) & meshR.Sf())
        + fvc::ddtPhiCorr(rUAR, UR, phiR)
        + phiQR;
    
    adjustPhi(phiR, UR, pR);
    
    // Non-orthogonal pressure corrector loop
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        // Pressure corrector
        fvScalarMatrix pEqnR(fvm::laplacian(rUAR, pR) == fvc::div(phiR));
        
        pEqnR.setReference(pRefCell, pRefValue);
        
        if (corr == nCorr-1 && nonOrth == nNonOrthCorr )
        {
            pEqnR.solve(meshR.solver("pRFinal"));
        }
        else
        {
            pEqnR.solve();
        }
        
        if (nonOrth == nNonOrthCorr)
        {
            phiR -= pEqnR.flux();
        }
    }
    
        #include "continuityErrsR.H"
    
    if(oCorr != nOuterCorr-1) { pR.relax(); }
    
    UR += fvc::reconstruct(phiQR) - rUAR*fvc::grad(pR);
    UR.correctBoundaryConditions();
}
