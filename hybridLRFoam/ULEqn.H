//      LES part.
fvVectorMatrix UEqnL
(
    fvm::ddt(UL)
    + fvm::div(phiL, UL)
    + sgsModel->divDevReff(UL)
    ==
    gradPL.flowDirection() * gradPL.value()
);

// alpha => upwinding factor
// alpha = 0 is central scheme; alpha = 1 is upwinding
if(blendConvection)
{
    if ( alphaPtr.valid() )
    {
        UEqnL += alphaPtr() * 
            ( 
                fvm::div(phiL, UL, "div(phiL,UL)|upwinding") 
                - fvm::div(phiL, UL)
            );
    } 
    else
    {
        FatalErrorIn("ULEqn.H ") 
            << "alphaPtr not valid! Check your blendConvection switch."
            << endl << abort(FatalError);
    }
}

if (oCorr == nOuterCorr-1)
{
    UEqnL.relax(1);
}
else
{
    UEqnL.relax();
}

volScalarField rUAL = 1.0/UEqnL.A();
surfaceScalarField ddtPhiQL = fvc::interpolate(QLes) & meshL.Sf();
surfaceScalarField phiQL = fvc::interpolate(rUAL) * ddtPhiQL;

if (momentumPredictor)
{
    solve(
        UEqnL == 
        ( -fvc::grad(pL) + fvc::reconstruct(ddtPhiQL) )
    );
}
