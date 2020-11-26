// p equation RANS part

for (int corr=0; corr<nCorr; corr++)
{
    U = rUA*UEqn.H();
    
    forAll(p.boundaryField(), patchi)
    {
        if (isA<zeroGradientFvPatchScalarField>(p.boundaryField()[patchi]))
        {
            phiQ.boundaryField()[patchi] = 0.0;
        }
    }
    
    phi = (fvc::interpolate(U) & mesh.Sf())
        + fvc::ddtPhiCorr(rUA, U, phi)
        + phiQ;
    
    adjustPhi(phi, U, p);
    
    // Non-orthogonal pressure corrector loop
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        // Pressure corrector
        fvScalarMatrix pEqn(fvm::laplacian(rUA, p) == fvc::div(phi));
        
        pEqn.setReference(pRefCell, pRefValue);
        
        if (corr == nCorr-1 && nonOrth == nNonOrthCorr )
        {
            pEqn.solve(mesh.solver("pRFinal"));
        }
        else
        {
            pEqn.solve();
        }
        
        if (nonOrth == nNonOrthCorr)
        {
            phi -= pEqn.flux();
        }
    }
    
        #include "continuityErrs.H"
    
    if(oCorr != nOuterCorr-1) { p.relax(); }
    
    U += fvc::reconstruct(phiQ) - rUA*fvc::grad(p);
    U.correctBoundaryConditions();
}
