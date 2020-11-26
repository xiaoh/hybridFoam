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

Application
    averageSpanwise

Description
    Post-processes data from flow simulation with homogeneous spanwise direction
    It look for "Mean" fields (generated via time averaging) and perform averaging on them.

    For each time: calculate: txx, txy,tyy, txy. Assuming that the
    mesh is periodic in the z directions, collapse Umeanx, Umeany, txx,
    txy and tyy to a plane and print them to fields. Automatic.

    Treatment of the boundary fields of the 2D fiels:
    Assume zeroGradient BC and extrapolate from internal fields.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "homogeneousSpanIndex.H"
#include "makeGraph.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "scalar.H"
#include "AverageFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    timeSelector::addOptions();

    #include "addRegionOption.H"

    argList::validOptions.insert("postfix", "postfix");  
    argList::validOptions.insert("fields", "\"(list of fields)\"");
    
    #include "setRootCase.H"
    #include "createTime.H"

    HashSet<word> selectedFields;
     if (args.optionFound("fields"))
     {
         args.optionLookup("fields")() >> selectedFields;
     }

     word postfix("");
     if (args.optionFound("postfix"))
     {
         args.optionLookup("postfix")() >> postfix;
     }


     // Get times list
     instantList timeDirs = timeSelector::select0(runTime, args);
    
     if (timeDirs.empty())
     {
         FatalErrorIn(args.executable())
             << "No times selected"
             << exit(FatalError);
     }
     
     #include "createNamedMesh.H"

     if (regionName != fvMesh::defaultRegion)
     {
         postfix = regionName;
     }

    #include "createMesh2D.H"

    Info << "Analyzing mesh regions ... " << endl;
    homogeneousSpanIndex homogeneousSpanIndexing(mesh);

    Info << "Checking mesh consistency ... " << nl << endl;
    vectorField CcMeanValues = homogeneousSpanIndexing.collapse(mesh.C());
    scalar Z = gAverage(mesh2D.C().component(vector::Z)());
    // Match fields
    forAll (CcMeanValues, i)
    {
        CcMeanValues[i].component(vector::Z) =  Z;
        label cellI = mesh2D.findNearestCell(CcMeanValues[i]);
        scalar dist = mag( mesh2D.C()[cellI] - CcMeanValues[i]);
        if (dist*dist*dist > 1e-12 * max(mesh.V()).value())
        {
            
            Info << "Z: " << Z << nl
                 << "Min distance: " << dist << nl 
                 << "Cell ID in 2D mesh: " << cellI << endl;
            FatalErrorIn("averageSpanwise")
                << "Cell " << i << " does not match any cell in the 2D field." << abort(FatalError);
        }
    }

    // For each time step read all fields
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Collapsing fields for time " << runTime.timeName() << endl;
        
        AverageFields<scalar>(mesh, mesh2D, homogeneousSpanIndexing, selectedFields);
        AverageFields<vector>(mesh, mesh2D, homogeneousSpanIndexing, selectedFields);
        AverageFields<symmTensor>(mesh, mesh2D, homogeneousSpanIndexing, selectedFields);
        AverageFields<tensor>(mesh, mesh2D, homogeneousSpanIndexing, selectedFields);
    }

    Info<< "\nEnd" << endl;

    return 0;
}


// ************************************************************************* //
