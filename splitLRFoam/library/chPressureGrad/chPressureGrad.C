/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
  -------------------------------------------------------------------------------
  License
  This file is part of OpenFOAM.

  OpenFOAM is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with OpenFOAM; if not, write to the Free Software Foundation,
  Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

  Class
  chPressureGrad

  \*---------------------------------------------------------------------------*/

#include "chPressureGrad.H"

namespace Foam {
  bool chPressureGrad::globalParamSet_ = false;
  Switch chPressureGrad::chFlowMode_ = false;
  vector chPressureGrad::flowDirection_ = vector::zero;
  dimensionedVector  chPressureGrad::Ubar_("Ubar", dimVelocity, vector::zero);
  dimensionedScalar  chPressureGrad::magUbar_ = mag(Ubar_);

  void  chPressureGrad::initPressureGrad(const IOdictionary& transportProperties)
  {
    
    if(transportProperties.found("Ubar"))
      {
        //  Read centerline velocity for channel simulations
        Ubar_ = transportProperties.lookup("Ubar");
        magUbar_ = mag(Ubar_);
        flowDirection_ = (Ubar_/magUbar_).value();
        chFlowMode_ = true;
        
        Info << "*********** NOTE ***********" << nl
             << "Running in CHANNEL FLOW MODE. " << nl
             << "Constant mass flux enforced!" << nl
             << "Ubar = " << Ubar_.value() << nl 
             << "****************************" << nl
             << endl;
      }
    else
      {
        Info << "+++++++++++ NOTE +++++++++++++"<< nl
             << "Running in REGULAR MODE."      << nl
             << "NO PRESSURE GRADIENT APPLIED!" << nl
             << "++++++++++++++++++++++++++++++"
             << endl;
      }
   
    globalParamSet_ = true;
  }

  // * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
  // * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

  // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
  chPressureGrad::chPressureGrad
  (
   volVectorField& U,
   word name,
   word solverName,
   word dictDir
   )
  :
    U_(U), 
    name_(name),
    solverName_(solverName),
    value_(name,  dimensionSet(0, 1, -2, 0, 0), 0.0),
    gradPDict_
    (
     IOobject
     (
      name,
      U.time().timeName(),
      dictDir,
      U.mesh(),
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
      )
     )
  {
    if(! globalParamSet_)
      FatalErrorIn("chPressureGrad::chPressureGrad()")
        << "Must first set global parameters by calling " 
        << "chPressureGrad::initPressureGrad(transportProperties)" << endl
        << abort(FatalError);

    if(chFlowMode_)
      {
        if(!gradPDict_.headerOk())
          {
            gradPDict_.add(name, value_);
          }
        else
          {
            gradPDict_.readIfPresent(name, value_);
          }
        
        Info << "Initializing " << name_ << " (from file or by default) as: "
             << name << " = " << value_.value() << endl;
      }
    else
      {
        gradPDict_.writeOpt() = IOobject::NO_WRITE; 
      }

  }
  
  
  // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

  // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

  void chPressureGrad::adjust(const volScalarField & rUA)
  {
    if(chFlowMode_)
      {
        dimensionedScalar magUbarStar =
          (flowDirection_ & U_)().weightedAverage(U_.mesh().V());
        
        // Calculate the pressure gradient increment needed to
        // adjust the average flow-rate to the correct value;
        // magUbar: target; magUbarStar: actual in simulation.
        dimensionedScalar gragPplus =
          (magUbar_ - magUbarStar)/rUA.weightedAverage(U_.mesh().V());
        
        U_ += flowDirection_ * rUA * gragPplus;
        value_ += gragPplus;
        gradPDict_.set(name_, value_);
        
        Info<< solverName_ 
            << " uncorrected Ubar = " << magUbarStar.value() 
            << "  "
            << name_ << " = " << value_.value() << endl;
      }
  }

  // * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


  // * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


  // * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

} // namespace Foam


// ************************************************************************* //
