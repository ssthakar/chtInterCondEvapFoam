if (finalIter)
{
    mesh.data::add("finalIteration", true);
}
// we don't really need this as fluid region and solid region are closely coupled
//  change in temp field in solid region does effect the fluid region
if (frozenFlow)
{
    #include "TEqn.H"

}


// correct for phase change and other stuff in the first pimple outer loop
if(oCorr == 0)
{
  volScalarField divU("divU", fvc::div(fvc::absolute(phi, U)));
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
    MRF.update();//added
    if (correctPhi)
    {
        // Calculate absolute flux
        // from the mapped surface velocity
        phi = mesh.Sf() & Uf();
        #include "correctPhi.H"
          
        // Make the flux relative to the mesh motion
        fvc::makeRelative(phi, U);
    }
    if(checkMeshCourantNo)
    {
      #include "meshCourantNo.H" // added
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

mixture.correct();
#include "alphaControls.H" //add actual file instead of linking vof??
#include "alphaEqnSubCycle.H" // same as above??1
interface.correct();

#include "UEqn.H"
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
    rho = alpha1*rho1+alpha2*rho2; //update rho

}



