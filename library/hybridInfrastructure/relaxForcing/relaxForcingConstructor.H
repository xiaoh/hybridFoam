  // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
  relaxForcing::relaxForcing
  (  
      const volVectorField& UL,
      const volVectorField& UR,
      Foam::incompressible::muLESModel& sgsModelRef,
      Foam::incompressible::muRASModel& ransModelRef,
      volVectorField& UAvg,
      volSymmTensorField& RAvg,
      volSymmTensorField& RResolvedAvg,
      volSymmTensorField& RSgsAvg,
      volScalarField& epsilonAvg,
      volVectorField& QLes,
      volVectorField& QRans,
      volScalarField& mrEpsilonAvg,
      volSymmTensorField& mrRAvg,
      IOdictionary& relaxParameters
  )
    :
      UL_(UL), UR_(UR),
      sgsModelRef_(sgsModelRef), ransModelRef_(ransModelRef),
      UAvg_(UAvg), RAvg_(RAvg), 
      RResolvedAvg_(RResolvedAvg), RSgsAvg_(RSgsAvg),
      epsilonAvg_(epsilonAvg),

      QLes_(QLes), QRans_(QRans),

      runTime_(UL.time()),
      meshL_(UL.mesh()),
      meshR_(UR.mesh()), 
      nCellsGL_(returnReduce(meshL_.nCells(), sumOp<label>())),
      nCellsGR_(returnReduce(meshR_.nCells(), sumOp<label>())),

      // Dictionaries
      relaxParameters_(relaxParameters),
      timeScaleDict_(relaxParameters.subDictPtr("timeScales")),
      couplingDict_(relaxParameters.subDictPtr("couplingOptions")),
      
      QSgs_
      (
          IOobject
          (
              "QSgs",
              runTime_.timeName(),
              meshL_,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE
          ),
          QLes_ * 0.0
      ),

      QSgsAvg_
      (
          IOobject
          (
              "QSgsAvg",
              runTime_.timeName(),
              meshL_,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE
          ),
          QSgs_
      ),

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

      MRelaxTime_
      (
          dimensioned<scalar>::lookupOrAddToDict
          ("MRelaxTime", timeScaleDict_, 
          averagingTime_.value() * scalar(0.2),
          dimTime)
      ),

      turbRelaxTime_
      (
          dimensioned<scalar>::lookupOrAddToDict
          ("turbRelaxTime", timeScaleDict_, 
          averagingTime_.value() * scalar(0.05),
          dimTime)
      ),

      nu0_(ransModelRef_.transport().lookup("nu")),
      nuRelaxFactor_(  
          dimensioned<scalar>::lookupOrAddToDict
          ("nuRelaxFactor", timeScaleDict_, 
          scalar(3.0),
          dimless)
      ),
      nuRelaxFactorL_(  
          dimensioned<scalar>::lookupOrAddToDict
          ("nuRelaxFactorL", timeScaleDict_, 
          nuRelaxFactor_.value(),
          dimless)
      ),
      nuRelaxFactorR_(  
          dimensioned<scalar>::lookupOrAddToDict
          ("nuRelaxFactorR", timeScaleDict_, 
          nuRelaxFactor_.value(),
          dimless)
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
      flexMask_
      (
          couplingDict_.lookupOrAddDefault<symmTensor>
          (
              "flexMask",
              symmTensor(0, 1, 0, 0, 0, 0)
          )
      ),
      epsilonSgsBy_
      (
          couplingDict_.lookupOrAddDefault<word>
          (
              "epsilonSgsBy", "dudx" 
          )
      ),
      rmDivTauSgs_
      (
          couplingDict_.lookupOrAddDefault<Switch>
          (
              "rmDivTauSgs", false
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
              meshL_,
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
      mapR2LEvery_
      (
          couplingDict_.lookupOrAddDefault<label>("mapR2LEvery", 1)
      ),
      mapL2REvery_
      (
          couplingDict_.lookupOrAddDefault<label>("mapL2REvery", 1)
      ),
      
      couplingStatus_(false),
      neutralQMomentum_(couplingDict_.lookupOrAddDefault<Switch>("neutralQMomentum", false)),

      ksmall_("ksmall", sqr(dimVelocity), 10*SMALL),
      vsmall_("vsmall", dimVolume, 10*SMALL),
      kSgsAvg_
      (
          IOobject
          (
              "kSgsAvg",
              runTime_.timeName(),
              meshL_,
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
              meshL_,
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
              meshL_,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE
          ),
          sgsModelRef_.k()
      ),
      RCorrPtr_(NULL),

      turbFlag_(resolutionModel::New(meshL_, meshR_)),
      lesFlagL_(turbFlag_ -> lesFlagL()),
      lesFlagR_(turbFlag_ -> lesFlagR()),

      meshTalk_ (meshL_, meshR_, mrEpsilonAvg, mrRAvg),

      consistencyDiagnosis_(couplingDict_.lookupOrAddDefault<Switch>("consistencyDiagnosis", false)),
      forcingDiagnosis_(couplingDict_.lookupOrAddDefault<Switch>("forcingDiagnosis", false))
{
    // Bound the coupling time parameters
    if ( couplingStartTime_ < runTime_.deltaT() || gradualCouplingDuration_ < runTime_.deltaT() )
    {
        FatalErrorIn("relaxForcing constructor: ") 
            << "Improperly specified time for coupling: "
            << "couplingStartTime = " << couplingStartTime_.value() << nl
            << "gradualCouplingDuation = " << gradualCouplingDuration_.value()
            << endl << abort(FatalError);
    }

    wordRe wrRConsisLevel(RConsisLevel_);
    wordRe wrPattern("[A-Za-z0-9_]*Direct", wordRe::DETECT);
    if (wrPattern.match(wrRConsisLevel))
    {
        Info << "Direct forcing strategy specified. I will compute corrective stresses." << endl;
        RCorrPtr_.reset(
            new volSymmTensorField
            (
                IOobject
                (
                    "RCorr",
                    runTime_.timeName(),
                    meshL_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                meshL_,
                dimensionedSymmTensor
                ("RCorr", dimVelocity*dimVelocity, symmTensor::zero)
            )
        );
    }
    else
    {
        Info << "Consisent forcing specified. No corrective stresses." << endl;
    }

    // Disable writeing QSgs and QSgsAvg unless used with impactFoam
    if ( rmDivTauSgs_ ) 
    {
        Info << "divTauSgs will be subtracted from QLes. "
             << "Are you running impactFoam with ADMRT sgs model? I am checking ..." << endl;
        if(epsilonSgsBy_ == "kdelta")
        {
            Info << "Yes. Then you are fine!" << endl;
        }
        else
        {
            FatalErrorIn("relaxForcing::relaxForcing()")
                << "Oops! Inconsistency found! " << nl
                << "Turn on option rmDivTauSgs ONLY IF you are running impactFoam with ADMRT model."
                << abort(FatalError);
        }
    }
    else
    {
        // QSgs_.writeOpt() = IOobject::NO_WRITE;
        // QSgsAvg_.writeOpt() = IOobject::NO_WRITE;
        QSgs_ *= 0.0;
        QSgsAvg_ *= 0.0;
    }

    checkRestart();
    initConsistentDict();
    printRelaxParameters();
}

