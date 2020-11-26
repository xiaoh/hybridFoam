    scalarField UMeanXvalues = channelIndexing.collapse
    (
        UMean.component(vector::X)()
    );

    scalarField UMeanYvalues = channelIndexing.collapse
    (
        UMean.component(vector::Y)()
    );

    scalarField UMeanZvalues = channelIndexing.collapse
    (
        UMean.component(vector::Z)()
    );

    scalarField RxxValues = channelIndexing.collapse(Rxx);
    scalarField RyyValues = channelIndexing.collapse(Ryy);
    scalarField RzzValues = channelIndexing.collapse(Rzz);
    scalarField RxyValues = channelIndexing.collapse(Rxy, true);

    scalarField pPrime2MeanValues = channelIndexing.collapse(pPrime2Mean);

    scalarField urmsValues = sqrt(mag(RxxValues));
    scalarField vrmsValues = sqrt(mag(RyyValues));
    scalarField wrmsValues = sqrt(mag(RzzValues));

    scalarField kValues =
        0.5*(sqr(urmsValues) + sqr(vrmsValues) + sqr(wrmsValues));

    const scalarField& y = channelIndexing.y();

    makeGraph(y, UMeanXvalues, "Uf", UMean.path(), gFormat);
    makeGraph(y, urmsValues, "u", UMean.path(), gFormat);
    makeGraph(y, vrmsValues, "v", UMean.path(), gFormat);
    makeGraph(y, wrmsValues, "w", UMean.path(), gFormat);
    makeGraph(y, RxyValues, "uv", UMean.path(), gFormat);
    makeGraph(y, kValues, "k", UMean.path(), gFormat);

    makeGraph(y, pPrime2MeanValues, "pPrime2Mean", UMean.path(), gFormat);

    // Compute and plot SGS contributions
   if (BMean.headerOk())
    {
      // SGS stress components
      scalarField BxxValues = channelIndexing.collapse(Bxx);
      scalarField ByyValues = channelIndexing.collapse(Byy);
      scalarField BzzValues = channelIndexing.collapse(Bzz);
      scalarField BxyValues = channelIndexing.collapse(Bxy, true);
      
      // Sgs rms velocities
      scalarField urmsSgsValues = sqrt(mag(BxxValues));
      scalarField vrmsSgsValues = sqrt(mag(ByyValues));
      scalarField wrmsSgsValues = sqrt(mag(BzzValues));
      
      scalarField kSgsValues =
        0.5*(sqr(urmsSgsValues) + sqr(vrmsSgsValues) + sqr(wrmsSgsValues));
      
      makeGraph(y, urmsSgsValues, "uB", UMean.path(), gFormat);
      makeGraph(y, vrmsSgsValues, "vB", UMean.path(), gFormat);
      makeGraph(y, wrmsSgsValues, "wB", UMean.path(), gFormat);
      makeGraph(y, BxyValues, "uvB", UMean.path(), gFormat);
      makeGraph(y, kSgsValues, "kB", UMean.path(), gFormat);
    }
    else
    {
      Info << "No BMean field found!" << endl;
    }
