//          |              [40 spaces]             |
//          v (first character)                    v (last space)
 Info << "\nhybridLRFoam version:                   " << hybridLRbuild << endl;;
relaxForcing::printVersion();
resolutionModel::printVersion();
Info << endl;

// Steps towards version printing (for executable; using hybridLRFoam as example)
// 1. copy and adapt hybridLRFoamVerson.{H,Cver} accordingly
// 2. copy and adapt printVersions.H 
// 3. Modify Make/file and Make/options
// 4. Add "#include printVersions.H" and "hybridLRFoamVerson.H" in hybridLRFoam.C
