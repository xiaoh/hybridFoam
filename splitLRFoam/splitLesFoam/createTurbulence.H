// create turbulence models and relaxation forcing

autoPtr<incompressible::muLESModel> sgsModel
(
    incompressible::muLESModel::New
    (
        U, phi, laminarTransport, RAvg, epsilonAvg
    )
);


// Create relaxForcing instance
lesForcing relaxQ
(
    U,
    sgsModel(),  // Extract reference from autoPtr type
    UAvg, RAvg, RResolvedAvg, RSgsAvg, epsilonAvg, 
    QLes, relaxParameters
);

