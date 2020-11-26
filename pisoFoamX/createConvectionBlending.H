IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);


dimensionedScalar nu
(
    transportProperties.lookup("nu")
);

autoPtr<volScalarField> alphaPtr(NULL);

Switch blendConvection
(
    transportProperties.lookupOrAddDefault<Switch>("blendConvection", false)
);

if(blendConvection)
{ 
    Info <<"Creating alpha ..." << endl;

    alphaPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "alpha",
                    runTime.constant(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("azero", dimless, 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );

    boundBox vortexInside
        (
            transportProperties.lookup("vortexInside")
        );

    vector blendingRange
        (
            transportProperties.lookup("blendingRange")
        );

    boundBox innerBox (vortexInside);
    boundBox outerBox (vortexInside);

    outerBox.min() -= blendingRange;
    outerBox.max() += blendingRange;
    
    DynamicList<label> outerCells;
    DynamicList<label> innerCells;
    
    const  volVectorField & cellCenters = mesh.C();
    forAll(cellCenters, celli)
    {
        point pt = cellCenters[celli];
        if (! outerBox.contains(pt)) // outside the outerBox
        {
            outerCells.append(celli);
        }
        else if ( innerBox.contains(pt) ) // inside the innerBox
        {
            innerCells.append(celli);
        }
    }
    
    labelList centralSchemeCells(innerCells.shrink().xfer());
    labelList upwindingCells(outerCells.shrink().xfer());
    
    fvScalarMatrix alphaEqn( fvm::laplacian(alphaPtr()) );
    scalarField onesTmpField(upwindingCells.size(), 1.0);
    scalarField zerosTmpField(centralSchemeCells.size(), 0.0);
    alphaEqn.setValues(upwindingCells, onesTmpField);
    alphaEqn.setValues(centralSchemeCells, zerosTmpField);

    solve(alphaEqn);
    alphaPtr().write();
}
