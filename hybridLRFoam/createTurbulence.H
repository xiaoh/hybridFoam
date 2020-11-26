// create turbulence models and relaxation forcing

autoPtr<incompressible::muLESModel> sgsModel
(
    incompressible::muLESModel::New
    (
        UL, phiL, laminarTransportL, RAvg, epsilonAvg
    )
);

autoPtr<incompressible::muRASModel> ransModel
(
    incompressible::muRASModel::New
    (
        UR, phiR, laminarTransportR, mrRAvg, mrEpsilonAvg
    )
);


autoPtr<volSymmTensorField> RPtr(NULL);
if(! meshR.foundObject<const volSymmTensorField>("RR"))
{
    RPtr.reset(
        new volSymmTensorField
        (
            IOobject
            (
                "RR",
                runTime.timeName(),
                meshR,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            meshR,
            dimensionedSymmTensor
            ("RR", dimVelocity*dimVelocity, symmTensor::zero)
        )
    );
}

// Create relaxForcing instance
relaxForcing relaxQ
(
    UL, UR,
    sgsModel(), ransModel(),  // Extract reference from autoPtr type
    UAvg, RAvg, RResolvedAvg, RSgsAvg, epsilonAvg, 
    QLes, QRans, mrEpsilonAvg, mrRAvg, relaxParameters
);

