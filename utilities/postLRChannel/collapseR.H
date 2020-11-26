    scalarField UMeanXvalues = channelIndexingR.collapse
    (
        UMean.component(vector::X)()
    );

    scalarField UMeanYvalues = channelIndexingR.collapse
    (
        UMean.component(vector::Y)()
    );

    scalarField UMeanZvalues = channelIndexingR.collapse
    (
        UMean.component(vector::Z)()
    );

    scalarField kMeanValues = channelIndexingR.collapse
    (
        kRMean
    );

/*
    scalarField RxxValues = channelIndexingR.collapse(Rxx);
    scalarField RyyValues = channelIndexingR.collapse(Ryy);
    scalarField RzzValues = channelIndexingR.collapse(Rzz);
    scalarField RxyValues = channelIndexingR.collapse(Rxy, true);

    scalarField pPrime2MeanValues = channelIndexingR.collapse(pPrime2Mean);


    scalarField epsilonValues = channelIndexingR.collapse(epsilonMean);

    scalarField nuMeanValues = channelIndexingR.collapse(nuMean);
    scalarField nuPrimeValues = channelIndexingR.collapse(nuPrime);

    scalarField gammaDotMeanValues = channelIndexingR.collapse(gammaDotMean);
    scalarField gammaDotPrimeValues = channelIndexingR.collapse(gammaDotPrime);


    scalarField urmsValues = sqrt(mag(RxxValues));
    scalarField vrmsValues = sqrt(mag(RyyValues));
    scalarField wrmsValues = sqrt(mag(RzzValues));
*/



    const scalarField& y = channelIndexingR.y();

    makeGraph(y, UMeanXvalues, "UfR", UMean.path(), gFormat);
    makeGraph(y, kMeanValues, "kR", UMean.path(), gFormat);

   if (RMean.headerOk())
    {
      scalarField RRxxValues = channelIndexingR.collapse(RRxx);
      scalarField RRyyValues = channelIndexingR.collapse(RRyy);
      scalarField RRzzValues = channelIndexingR.collapse(RRzz);
      scalarField RRxyValues = channelIndexingR.collapse(RRxy, true);
      
      makeGraph(y, RRxxValues, "Rxx", UMean.path(), gFormat);
      makeGraph(y, RRyyValues, "Ryy", UMean.path(), gFormat);
      makeGraph(y, RRzzValues, "Rzz", UMean.path(), gFormat);
      makeGraph(y, RRxyValues, "Rxy", UMean.path(), gFormat);
    }
    else
    {
      Info << "No RMean field found!" << endl;
    }


//    makeGraph(y, urmsValues, "uR", UMean.path(), gFormat);
//    makeGraph(y, vrmsValues, "vR", UMean.path(), gFormat);
//    makeGraph(y, wrmsValues, "wR", UMean.path(), gFormat);
//    makeGraph(y, RxyValues, "uvR", UMean.path(), gFormat);


//    makeGraph(y, pPrime2MeanValues, "pPrime2MeanR", UMean.path(), gFormat);

    /*
    makeGraph(y, epsilonValues, "epsilon", UMean.path(), gFormat);
    makeGraph(y, nuMeanValues, "nu", UMean.path(), gFormat);
    makeGraph(y, nuPrimeValues, "nuPrime", UMean.path(), gFormat);
    makeGraph(y, gammaDotMeanValues, "gammaDot", UMean.path(), gFormat);
    makeGraph(y, gammaDotPrimeValues, "gammaDotPrime", UMean.path(), gFormat);
    */
