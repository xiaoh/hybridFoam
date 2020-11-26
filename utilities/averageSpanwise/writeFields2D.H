      Info << "Writing 2D fields ...\n" << endl;

      U2.write();
      pPrime2Mean2.write();
      RResolved2.write();

   if (BMean.headerOk())
    {
      B2.write();
    }
