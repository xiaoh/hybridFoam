// create turbulence models and relaxation forcing

autoPtr<incompressible::muRASModel> ransModel
(
    incompressible::muRASModel::New
    (
        U, phi, laminarTransport, mrRAvg, mrEpsilonAvg
    )
);


autoPtr<volSymmTensorField> RPtr(NULL);
if(! mesh.foundObject<const volSymmTensorField>("RR"))
{
    RPtr.reset(
        new volSymmTensorField
        (
            IOobject
            (
                "RR",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedSymmTensor
            ("RR", dimVelocity*dimVelocity, symmTensor::zero)
        )
    );
}

// Create relaxForcing instance
ransForcing relaxQ
(
    U,
    ransModel(),  // Extract reference from autoPtr type
    QRans, mrEpsilonAvg, mrRAvg, relaxParameters
);

