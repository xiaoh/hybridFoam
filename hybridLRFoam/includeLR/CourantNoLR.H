/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

Global
    CourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

\*---------------------------------------------------------------------------*/

scalar CoNumLES = 0.0;
scalar meanCoNumLES = 0.0;

if(solveLES)
{
    if (meshL.nInternalFaces())
    {
        surfaceScalarField SfUfbyDelta = 
            meshL.surfaceInterpolation::deltaCoeffs()*mag(phiL);
        
        CoNumLES = max(SfUfbyDelta/meshL.magSf())
            .value()*runTime.deltaT().value();
        
        meanCoNumLES = (sum(SfUfbyDelta)/sum(meshL.magSf()))
            .value()*runTime.deltaT().value();
    }
    
    Info<< "LES Courant Number mean: " << meanCoNumLES
        << " max: " << CoNumLES << endl;
}

scalar CoNumRans = 0.0;
scalar meanCoNumRans = 0.0;

if(solveRANS)
{
    if (meshR.nInternalFaces())
    {
        surfaceScalarField SfUfbyDelta = 
            meshR.surfaceInterpolation::deltaCoeffs()*mag(phiR);
        
        CoNumRans = max(SfUfbyDelta/meshR.magSf())
            .value()*runTime.deltaT().value();
        
        meanCoNumRans = (sum(SfUfbyDelta)/sum(meshR.magSf()))
            .value()*runTime.deltaT().value();
    }
    
    Info<< "RAN Courant Number mean: " << meanCoNumRans
        << " max: " << CoNumRans << endl;
}
Info << "---" << endl;
// ************************************************************************* //
