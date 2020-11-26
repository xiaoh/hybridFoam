//
// createMesh2D.H
// ~~~~~~~~~~~~

    Foam::Info
        << "Create 2D mesh for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;

word Case2D("Case2D"+postfix);

if (! isDir(Case2D))
  {
      Info << Case2D << " not found! You should create it using foamTwoD" << endl;
  }

Foam::fvMesh mesh2D
(
    Foam::IOobject
    (
        Case2D,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);

