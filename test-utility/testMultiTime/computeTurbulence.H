// Solver turbulence equations

 sgsModel->correct();
 ransModel->correct();

 relaxQ.correct();

 // update B (only for the averaging, not part of actual computing)
 B = sgsModel->B();

   
