/*---------------------------------------------------------------------------* \
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

Application
    Lambda2

Description
    Calculates and writes the second largest eigenvalue of the sum of the
    square of the symmetrical and anti-symmetrical parts of the velocity
    gradient tensor.

    The -noWrite option has no meaning.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    IOobject Uheader
    (
        "UL",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (Uheader.headerOk())
    {
        Info<< "    Reading UL" << endl;
        volVectorField U(Uheader, mesh);

        const volTensorField gradU(fvc::grad(U));

        volTensorField SSplusWW =
            (symm(gradU) & symm(gradU)) + (skew(gradU) & skew(gradU));

        volScalarField Lambda2
        (
            IOobject
            (
                "Lambda2",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            -eigenValues(SSplusWW)().component(vector::Y)
        );

        Info << "    Writing -Lambda2" << endl;
        Lambda2.write();
    }
    else
    {
        Info<< "    Looked for UL but not found" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
