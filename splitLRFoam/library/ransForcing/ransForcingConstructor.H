  // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
  ransForcing::ransForcing
  (  
      const volVectorField& U,
      Foam::incompressible::muRASModel& ransModelRef,
      volVectorField& QRans,
      volScalarField& mrEpsilonAvg,
      volSymmTensorField& mrRAvg,
      IOdictionary& relaxParameters
  )
    :
      U_(U),
      ransModelRef_(ransModelRef),
      QRans_(QRans),
      
      runTime_(U.time()),
      mesh_(U.mesh()), 
      nCellsG_(returnReduce(mesh_.nCells(), sumOp<label>())),

      // Dictionaries
      relaxParameters_(relaxParameters),
      timeScaleDict_(relaxParameters.subDictPtr("timeScales")),
      couplingDict_(relaxParameters.subDictPtr("couplingOptions")),
      
      // Time scales:
      MRelaxTime_
      (
          dimensioned<scalar>::lookupOrAddToDict
          ("MRelaxTime", timeScaleDict_, 
          1.0,
          dimTime)
      ),

      turbRelaxTime_
      (
          dimensioned<scalar>::lookupOrAddToDict
          ("turbRelaxTime", timeScaleDict_, 
          1.0,
          dimTime)
      ),

    // Coupling options:
      LesRansCoupling_
      (
          couplingDict_.lookupOrAddDefault<Switch>("LesRansCoupling", true)
      ),

      couplingStartTime_
      (
          dimensioned<scalar>::lookupOrAddToDict
          ("couplingStartTime", couplingDict_, 
          runTime_.endTime().value() * scalar(0.1),
          dimTime)
      ),

      gradualCouplingDuration_
      (
          dimensioned<scalar>::lookupOrAddToDict
          ("gradualCouplingDuration", couplingDict_,
          couplingStartTime_.value(), 
          dimTime)
      ),

      consistentDict_
      (
          IOobject
          (
              "consistentProperties",
              runTime_.timeName(),
              mesh_,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE,
              true
          )
      ),

      rampQ_(0),
      enforceMode_
      (
          couplingDict_.lookupOrAddDefault<word>
          (
              "enforceMode", "both" 
          )
      ),

      // Loose coupling dict entries
      mapL2REvery_
      (
          couplingDict_.lookupOrAddDefault<label>("mapL2REvery", 1)
      ),
      
      cleanFootprint_
      (
          couplingDict_.lookupOrAddDefault<Switch>("cleanFootprint", true)
      ),

      couplingStatus_(false),

      ksmall_("ksmall", sqr(dimVelocity), 10*SMALL),
      vsmall_("vsmall", dimVolume, 10*SMALL),

      lesFlag_
      (
          IOobject
          (
              "lesFlagR",
              runTime_.timeName(),
              mesh_,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE
          ),
          mesh_,
          dimensionedScalar("lfvalue", dimless, 0)
      ),

      ransRange_
        (
            "ransRange", 
            dimLength,
            couplingDict_.lookup("ransRange")
        ),

      mrUAvg_
      (
          IOobject
          (
              "mrUAvg",
              runTime_.timeName(),
              mesh_,
              IOobject::READ_IF_PRESENT,
              IOobject::NO_WRITE
          ),
          mesh_,
          dimensionedVector("mruavg", dimVelocity, vector(0, 0.00035, 0))
      ),
      mrEpsilonAvg_(mrEpsilonAvg),
      mrRAvg_(mrRAvg)

{
    // Bound the coupling time parameters
    if ( couplingStartTime_ < runTime_.deltaT() || gradualCouplingDuration_ < runTime_.deltaT() )
    {
        FatalErrorIn("ransForcing constructor: ") 
            << "Improperly specified time for coupling: "
            << "couplingStartTime = " << couplingStartTime_.value() << nl
            << "gradualCouplingDuation = " << gradualCouplingDuration_.value()
            << endl << abort(FatalError);
    }

    wallDist y(mesh_, true);
    lesFlag_.internalField() = pos(y.y() - ransRange_);

    checkRestart();
    initConsistentDict();
    printRelaxParameters();
}

