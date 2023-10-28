if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

if (frozenFlow)
{
    #include "TEqn.H"
}
// correct for phase change and other stuff in the first pimple outer loop
if(oCorr == 0)
{
  if(debug)
  {
    Info << " correcting for phase change and MULES"<<endl;
  }
  fluid.correctMassSources(T);
  fluid.solve();
  rho = fluid.rho();
}

if(debug)
{
  Info<< "printing out sizes of fields before mesh update\n"<<endl;
  Info<<"rho: "<<rho.size()<<endl;
  Info<<"phi: "<<phi.size()<<endl;
  Info<<"U: "<<U.size()<<endl;
  Info<<"rhoPhi: "<<rhoPhi.size()<<endl;
}

// adaptive mesh refinement and face flux update
if(oCorr == 0 || moveMeshOuterCorrectors)
{ 

  Info << "\nupdating mesh \n"<<endl;
  mesh.update();
  if(mesh.changing())
  {
    gh = (g & mesh.C()) - ghRef; 
    ghf = (g & mesh.Cf()) - ghRef;

    if (correctPhi)
    {
        // Calculate absolute flux
        // from the mapped surface velocity
        phi = mesh.Sf() & Uf();
        #include "correctPhi.H"
          
        // Make the flux relative to the mesh motion
        fvc::makeRelative(phi, U);
    }
  }
}

if(debug)
{
  Info<< "printing out sizes of fields after mesh update\n"<<endl;
  Info<<"rho: "<<rho.size()<<endl;
  Info<<"phi: "<<phi.size()<<endl;
  Info<<"U: "<<U.size()<<endl;
  Info<<"rhoPhi: "<<rhoPhi.size()<<endl;
}

if(debug)
{
  Info << "governing Equations being solved"<<endl;
}
  // solve the Equations
#include "UEqn.H"
#include "YEqns.H"
#include "TEqn.H" 
//- pressure corrector loop
for(int corr=0;corr<nCorr;corr++)
{
    if(debug)
    {
      Info << "pressure Equation called"<<endl;
    }
    #include "pEqn.H"
}
// correct turbulence
turbulence.correct(); 
// stop loop if finalIter
if (finalIter)
{
    mesh.data::remove("finalIteration");
    rho = fluid.rho();

}



