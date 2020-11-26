    // Last update:
    string versionLRFoam("Time-stamp: <2010-03-30 15:52:28 xiaoh>");

    Info<< "LES: Reading field p\n" << endl;
    volScalarField pL
    (
        IOobject
        (
            "pL",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "RAS: Reading field p\n" << endl;
    volScalarField pR
    (
        IOobject
        (
            "pR",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pL
    );


    Info<< "LES: Reading field U\n" << endl;
    volVectorField UL
    (
        IOobject
        (
            "UL",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

   
    Info<< "RAS: Reading/creating field U\n" << endl;
    volVectorField UR
    (
        IOobject
        (
            "UR",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        UL
    );

  volSymmTensorField B
    (
     IOobject
     (
      "B",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh,
     dimensionedSymmTensor
       ("B", dimVelocity*dimVelocity, symmTensor::zero)
     );

   
#   include "createPhiLR.H"

    label pRefCell = 0;
    scalar pRefValue = 0.0;

   // Read in pLRefcell & pLRefValue (also for R) and set 
   // the value pRefCell and pRefValue here.
   // The first setRefCell is dummy (ensure consistency)

   // setRefCell(pR, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
    setRefCell(pL, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

    singlePhaseTransportModel laminarTransportL(UL, phiL);
    singlePhaseTransportModel laminarTransportR(UR, phiR);


#   include "createRelaxFields.H"

    autoPtr<incompressible::LESModel> sgsModel
    (
     incompressible::LESModel::New(UL, phiL, laminarTransportL)
    );

    
    volScalarField* lesFlagPtr=NULL;
    
   //- Potential inconsistency between here and relaxForcing class.
   //- Need to be done in a more elegant way
   //- For now, the code need to be updated consistently anytime 
   //  wheevern a new resolution criteron is introduced.
    {
      // resolutionCriterion valid only in the scope of {}

      dictionary couplingDict(relaxParameters.subDictPtr("couplingOptions"));

      word resolutionCriterion
      (
       couplingDict.lookupOrAddDefault<word>
       (
        "resolutionCriterion", "length" 
        )
       );
    
      if(
       resolutionCriterion == "energy" || 
       resolutionCriterion == "Energy")
      {
        lesFlagPtr = & lesFlagK;
      }
    else if 
      (
       resolutionCriterion == "length" || 
       resolutionCriterion == "Length"
       )
      {
        lesFlagPtr = & lesFlagL;
      }
    else if
      (
       resolutionCriterion == "wallDist" || 
       resolutionCriterion == "walldist"
       )
      {
        lesFlagPtr = & lesFlagGeneral;
      }
    else
      {
        Info << "***Warning: createFieldsLR.H:: unknown resolutionCriterion "
             << "in dictionary relaxationParameters.\n "
             << "Possible values: Length/Energy/wallDist"
             << "Default to Length criterion." << endl;        

        lesFlagPtr = & lesFlagL;
      }

      if (! lesFlagPtr)
        FatalErrorIn("createFieldsLR.H")
          <<"lesFlagPtr not properly set."<< abort(FatalError);
      
      Info << "createFieldsLR.H:: Resolution criterion used: "
           << lesFlagPtr->name() << "\n" << endl;
    }


    autoPtr<incompressible::muRASModel> ransModel
    (
     incompressible::muRASModel::New(UR, phiR, laminarTransportR,
                                     RAvg, epsilonAvg, *lesFlagPtr)
    );

// Create relaxForcing instance
relaxForcing relaxQ
(
 UL, UR,
 sgsModel(), ransModel(),  // Extract reference from autoPtr type
 UAvg, RAvg, RResolvedAvg, RSgsAvg, epsilonAvg, turbLength, 
 lesFlagL,  lesFlagK, lesFlagGeneral, resIndicatorL, resIndicatorK,
 QLes, QRans, relaxParameters
 );
