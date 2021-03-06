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

File:
    createPhiLR

Description
    Creates and initialises the relative face-flux field phi for both fields (Les & Rans).

\*---------------------------------------------------------------------------*/

#ifndef createPhiLR_H
#define createPhiLR_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "Reading/calculating face flux field phi for LES\n" << endl;

surfaceScalarField phiL
(
    IOobject
    (
        "phiL",
        runTime.timeName(),
        meshL,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(UL) & meshL.Sf()
);

Info<< "RAS: Reading/calculating face flux field phi\n" << endl;

surfaceScalarField phiR
(
    IOobject
    (
        "phiR",
        runTime.timeName(),
        meshR,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(UR) & meshR.Sf()
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
