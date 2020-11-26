//      RANS part.
fvVectorMatrix UEqnR
(
    fvm::ddt(UR)
    + fvm::div(phiR, UR)
    + ransModel->divDevReff(UR)
    ==
    gradPR.flowDirection()*gradPR.value()
);

if (oCorr == nOuterCorr-1)
{
    UEqnR.relax(1);
}
else
{
    UEqnR.relax(); 
}

volScalarField rUAR = 1.0/UEqnR.A();
surfaceScalarField ddtPhiQR = fvc::interpolate(QRans) & meshR.Sf();
surfaceScalarField phiQR = fvc::interpolate(rUAR) * ddtPhiQR;

if (momentumPredictor)
{
    solve(
        UEqnR == 
        ( -fvc::grad(pR) + fvc::reconstruct(ddtPhiQR) )
    );
}
