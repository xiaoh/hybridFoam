//      LES part.
fvVectorMatrix UEqn
(
    fvm::ddt(U)
    + fvm::div(phi, U)
    + sgsModel->divDevReff(U)
    ==
    gradP.flowDirection() * gradP.value()
);

if (oCorr == nOuterCorr-1)
{
    UEqn.relax(1);
}
else
{
    UEqn.relax();
}

volScalarField rUA = 1.0/UEqn.A();
surfaceScalarField ddtPhiQ = fvc::interpolate(QLes) & mesh.Sf();
surfaceScalarField phiQ = fvc::interpolate(rUA) * ddtPhiQ;

if (momentumPredictor)
{
    solve(
        UEqn == 
        ( -fvc::grad(p) + fvc::reconstruct(ddtPhiQ) )
    );
}
