

    vectorField UMeanValues = homogeneousSpanIndexing.collapse(UMean);

    scalarField pPrime2MeanValues = homogeneousSpanIndexing.collapse(pPrime2Mean);
    symmTensorField UPrime2MeanValues =  homogeneousSpanIndexing.collapse(UPrime2Mean);
    symmTensorField BMeanValues =  homogeneousSpanIndexing.collapse(BMean);

   /* scalarField xMean = homogeneousSpanIndexing.collapse
    (
     mesh.C().component(vector::X)()
     );

     string  fileNameUa(UMean.path()/"Uall.xy");
     homogeneousSpanIndexing.graphField(fileNameUa, UMeanAll);
    */
    vectorField CcMeanValues = homogeneousSpanIndexing.collapse(mesh.C());

   scalar Z = gAverage(mesh2D.C().component(vector::Z)());

// Match fields
    forAll (UMeanValues, i)
    {
      CcMeanValues[i].component(vector::Z) =  Z;
      label cellI = mesh2D.findNearestCell(CcMeanValues[i]);
      scalar dist = mag( mesh2D.C()[cellI] - CcMeanValues[i]);
      if (dist*dist*dist > 1e-6 * max(mesh.V()).value())
        {

          Info << "Z: " << Z << endl;

          Info << "Min distance: " 
               << dist << "Cell ID in 2D mesh: " << cellI << endl;

          FatalErrorIn("collapse(..)")
            << "Cell " << i << " does not match any cell in the 2D field." << abort(FatalError);
        }
      else // Assign the fields
        {
          U2.internalField()[cellI] = UMeanValues[i]; 
          pPrime2Mean2.internalField()[cellI] = pPrime2MeanValues[i]; 
          RResolved2.internalField()[cellI] = UPrime2MeanValues[i]; 
          B2.internalField()[cellI] = BMeanValues[i]; 
        }
      // Info << i << ": Nearest cell: " << cellI << "; Dist = " << dist << endl;
    }

