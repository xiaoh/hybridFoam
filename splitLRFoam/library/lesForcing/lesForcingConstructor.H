  // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
  lesForcing::lesForcing
  (  
      const volVectorField& U,
      Foam::incompressible::muLESModel& sgsModelRef,
      volVectorField& UAvg,
      volSymmTensorField& RAvg,
      volSymmTensorField& RResolvedAvg,
      volSymmTensorField& RSgsAvg,
      volScalarField& epsilonAvg,
      volVectorField& QLes,
      IOdictionary& relaxParameters
  )
    :
      U_(U),
      sgsModelRef_(sgsModelRef),
      UAvg_(UAvg), RAvg_(RAvg), 
      RResolvedAvg_(RResolvedAvg), RSgsAvg_(RSgsAvg),
      epsilonAvg_(epsilonAvg),

      QLes_(QLes),
      
      runTime_(U.time()),
      mesh_(U.mesh()),
      nCellsG_(returnReduce(mesh_.nCells(), sumOp<label>())),


      // Dictionaries
      relaxParameters_(relaxParameters),
      timeScaleDict_(relaxParameters.subDictPtr("timeScales")),
      couplingDict_(relaxParameters.subDictPtr("couplingOptions")),
      
      // Time scales:
      averagingTime_
      (
          dimensioned<scalar>::lookupOrAddToDict
          ("averagingTime", timeScaleDict_, 
          runTime_.endTime().value() * scalar(0.1),
          dimTime)
      ),

      URelaxTime_
      (
          dimensioned<scalar>::lookupOrAddToDict
          ("URelaxTime", timeScaleDict_, 
          averagingTime_.value() * scalar(0.2),
          dimTime)
      ),

      RRelaxTime_
      (
          dimensioned<scalar>::lookupOrAddToDict
          ("RRelaxTime", timeScaleDict_, 
          averagingTime_.value() * scalar(0.3),
          dimTime)
      ),

    // Coupling options:
      LesRansCoupling_
      (
          couplingDict_.lookupOrAddDefault<Switch>("LesRansCoupling", true)
      ),
      RConsisLevel_
      (
          couplingDict_.lookupOrAddDefault<word>
          (
              "RConsisLevel", "k" 
          )
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

      // Loose coupling dict entries
      mapR2LEvery_
      (
          couplingDict_.lookupOrAddDefault<label>("mapR2LEvery", 1)
      ),

      cleanFootprint_
      (
          couplingDict_.lookupOrAddDefault<Switch>("cleanFootprint", true)
      ),
      
      couplingStatus_(false),

      ksmall_("ksmall", sqr(dimVelocity), 10*SMALL),
      vsmall_("vsmall", dimVolume, 10*SMALL),
      kSgsAvg_
      (
          IOobject
          (
              "kSgsAvg",
              runTime_.timeName(),
              mesh_,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE
          ),
          sgsModelRef.k()
      ),
      kAvg_
      (
          IOobject
          (
              "kAvg",
              runTime_.timeName(),
              mesh_,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE
          ),
          0.5*tr(RAvg_)
      ),
      kLES_
      (
          IOobject
          (
              "kLES",
              runTime_.timeName(),
              mesh_,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE
          ),
          sgsModelRef_.k()
      ),
      
      lesFlag_
      (
          IOobject
          (
              "lesFlag",
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

      mlUR_
      (
          IOobject
          (
              "mlUR",
              runTime_.timeName(),
              mesh_,
              IOobject::READ_IF_PRESENT,
              IOobject::NO_WRITE
          ),
          mesh_,
          dimensionedVector("mlURvalue", dimVelocity, vector(0, 0, 0))
      ),
      mlRR_
      (
          IOobject
          (
              "mlRR",
              runTime_.timeName(),
              mesh_,
              IOobject::READ_IF_PRESENT,
              IOobject::NO_WRITE
          ),
          mesh_,
          dimensionedSymmTensor("mlRRvalue", dimVelocity*dimVelocity, symmTensor::zero)
      ),
      mlkR_
      (
          IOobject
          (
              "mlkR",
              runTime_.timeName(),
              mesh_,
              IOobject::READ_IF_PRESENT,
              IOobject::NO_WRITE
          ),
          mesh_,
          ksmall_
      )
{
    // Bound the coupling time parameters
    if ( couplingStartTime_ < runTime_.deltaT() || gradualCouplingDuration_ < runTime_.deltaT() )
    {
        FatalErrorIn("lesForcing constructor: ") 
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

