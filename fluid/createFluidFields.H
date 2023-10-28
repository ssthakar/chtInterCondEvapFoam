// Things to do 
//
//
//
// set fluid region fields update


// p_rgh
PtrList<volScalarField> p_rghFluid(fluidRegions.size()); //done 
//- U
PtrList<volVectorField> UFluid(fluidRegions.size()); //done 
// createPhi from interCondFoam
PtrList<surfaceScalarField> phiFluid(fluidRegions.size()); //done 
// e based Thermo for fluidRegions
PtrList<twoPhaseMixtureEThermo> thermoFluid(fluidRegions.size()); //done 
// Twophasemixture for fluid regions
PtrList<temperaturePhaseChangeTwoPhaseMixture> mixtureFluid(fluidRegions.size()); //done 
//rho1 and rho2
PtrList<dimensionedScalar> rho1Fluid(fluidRegions.size());
PtrList<dimensionedScalar> rho2Fluid(fluidRegions.size());
 // density for fluid mixture 
PtrList<volScalarField> rhoFluid(fluidRegions.size());
// interface properties
PtrList<interfaceProperties> interfaceFluid(fluidRegions.size());
// turbulence model for fluid regions
PtrList<incompressible::turbulenceModel> turbFluid(fluidRegions.size());
//KappaEff for fluid regions 
PtrList<volScalarField> kappaEffFluid(fluidRegions.size()); // what do?? here 
//- gh 
PtrList<volScalarField> ghFluid(fluidRegions.size());
//- ghRef
PtrList<dimensionedScalar> ghRefFluid(fluidRegions.size());
//- ghFlux
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
// hRef
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
//- p
PtrList<volScalarField> pFluid(fluidRegions.size());
//- T
PtrList<volScalarField> TFluid(fluidRegions.size()); 
// turbulent Prandtl Number
PtrList<dimensionedScalar> PrtFluid(fluidRegions.size());
//- mass flux
PtrList<surfaceScalarField> rhoPhiFluid(fluidRegions.size());
//- field kinetic energy
PtrList<volScalarField> KFluid(fluidRegions.size());
//- rhoCp
PtrList<volScalarField> rhoCpFluid(fluidRegions.size());
//- correct mass flux
List<scalar> initialMassFluid(fluidRegions.size());
//- frozen field?
List<bool> frozenFlowFluid(fluidRegions.size(), false);
//- correctPhi for adaptive mesh refinement
List<bool> correctPhiFluid(fluidRegions.size());
//- taken from icoReacting
List<label> pRefCellFluid(fluidRegions.size());
//- taken from icoReacting
List<scalar> pRefValueFluid(fluidRegions.size());
//- construct universal gravity, common to all fluid regions
const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);
//-PIMPLE control
PtrList<pimpleControl> pimpleFluid(fluidRegions.size());
//- fvOptions and MRF
PtrList<fv::options> fvOptionsFluid(fluidRegions.size());
PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
//- total cont error
List<scalar> cumulativeContErrFluid(fluidRegions.size());
//alphaPhiFluid for createAlphaFluxes
PtrList<surfaceScalarField> alphaPhiFluid(fluidRegions.size());
// MULES compressed flux if needed
PtrList<surfaceScalarField> alphaPhiUnFluid(fluidRegions.size());

//- populate fluid field pointer lists
forAll(fluidRegions, i)
{
    
    pimpleFluid.set
    (
        i,
        new pimpleControl(fluidRegions[i])
    );

    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;
 
    Info<< "    Adding p_rghFluid\n" << endl;
    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<<"      adding to face flux field phi\n"<<endl;
    phiFluid.set
    (
      i,
      new surfaceScalarField
      (
        IOobject
        (
          "phi",
          runTime.timeName(),
          fluidRegions[i],
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
        ),
        fvc::flux(UFluid[i])
      )
    );


    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );
    
    // add thermo to pointer list
    Info<<"   adding to e based Thermo\n"<<endl;
    thermoFluid.set
    (
      i,
      new twoPhaseMixtureEThermo
      (
        UFluid[i],
        phiFluid[i]
      )
    );
    // set pointer list to new mixture
    Info<<"Creating temperaturePhaseChangeTwoPhaseMixture and adding it to fluidRegion\n"<<endl;
    mixtureFluid.set
    (
      i,
      temperaturePhaseChangeTwoPhaseMixture::New(thermoFluid[i], fluidRegions[i])
    );

    Info<<"     Adding rhoFluid\n" << endl;
    rhoFluid.set
    (
      i,
      new volScalarField
      (
        IOobject
        (
          "rho",
          runTime.timeName(),
          fluidRegions[i],
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
        ),
        thermoFluid[i].alpha1()*thermoFluid[i].rho1()
        +thermoFluid[i].alpha2()*thermoFluid[i].rho2()
      )
    );
    rhoFluid[i].oldTime();
    
    // set mass flux across faces
    rhoPhiFluid.set
    (
      i,
      new surfaceScalarField
      (
        IOobject
        (
          "rhoPhi",
          runTime.timeName(),
          fluidRegions[i],
          IOobject::NO_READ,
          IOobject::NO_WRITE
        ),
        fvc::interpolate(rhoFluid[i])*phiFluid[i]
      )
    );
    
    //set interface 
    interfaceFluid.set
    (
      i,
      new interfaceProperties
      (
        thermoFluid[i].alpha1(),
        UFluid[i],
        thermoFluid[i].transportPropertiesDict()
      )
    );


    Info<<"adding to turbulence PtrList\n"<<endl;
    turbFluid.set
    (
      i,
      incompressible::turbulenceModel::New
      (
        UFluid[i],
        phiFluid[i],
        thermoFluid[i]
      ).ptr()
    );
    turbFluid[i].validate();
    // Info<< "    Adding to TFluid\n" << endl;
    // TFluid.set
    // (
    //     i,
    //     new volScalarField
    //     (
    //         IOobject
    //         (
    //             "T",
    //             runTime.timeName(),
    //             fluidRegions[i],
    //             IOobject::MUST_READ,
    //             IOobject::AUTO_WRITE
    //         ),
    //         fluidRegions[i]
    //     )
    // );
      //--------------------------------------------------------------------------------
      /* this code does what the 3 include lines below do
      #include "readGravitationalAcceleration.H"
      #include "readhRef.H"
      #include "gh.H"

      */ 
    
    Info<< "    Adding hRefFluid\n" << endl;
    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("hRef", dimLength, Zero)
        )
    );
 
    Info<< "    Adding ghRef\n" << endl;
    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRefFluid[i]
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );
    //- set ghRef
    ghRefFluid.set(i, new dimensionedScalar(ghRef));
    
    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (g & fluidRegions[i].C()) - ghRef
        )
    );

    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (g & fluidRegions[i].Cf()) - ghRef
        )
    );
    //-------------------------------------------------------------------------
    

    // Info<< "    Adding to multiphaseSystemFluid\n" << endl;
    //   
    // phaseSystemFluid.set(i, multiphaseSystem::New(fluidRegions[i]).ptr());

    
    // if(debug)
    // {
    // Info<< "    updating p based on p_rgh\n" << endl;
    // }
    // pFluid.set
    // (
    //   i,
    //   new volScalarField
    //   (
    //     IOobject
    //     (
    //       "p",
    //       runTime.timeName(),
    //       fluidRegions[i],
    //       IOobject::NO_READ,
    //       IOobject::AUTO_WRITE
    //     ),
    //     p_rghFluid[i] + rhoFluid[i]*ghFluid[i]
    //   )
    // );
    //
    pRefCellFluid[i] = 0;
    pRefValueFluid[i] = 0.0;

    // if(debug)
    // {
    //   Info << " after p update"<<endl;
    // }

    // changed pFluid to thermoFluid.p()
    setRefCell
    (
        thermoFluid[i].p(), // don't build p directly like icoReacting
        p_rghFluid[i],
        pimpleFluid[i].dict(),
        pRefCellFluid[i],
        pRefValueFluid[i]
    );

    if (p_rghFluid[i].needReference())
    {
        thermoFluid[i].p() += dimensionedScalar
        (
            "p",
            thermoFluid[i].p().dimensions(),
            pRefValueFluid[i] - getRefCellValue(thermoFluid[i].p(), pRefCellFluid[i])
        );
        p_rghFluid[i] = thermoFluid[i].p() - rhoFluid[i]*ghFluid[i];
    }
    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());
    fluidRegions[i].setFluxRequired(thermoFluid[i].alpha1().name());

    fvOptionsFluid.set
    (
        i,
        new fv::options
        (
            fluidRegions[i]
        )
    );
    // set MRF zones for fluid regions, wont ever use them lol but just in case
    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );
    //set turbulent Prandtl number
    PrtFluid.set
    (
      i,
      new dimensionedScalar
      (
        "Prt",
        dimless,
        thermoFluid[i].transportPropertiesDict()
      )
    );
    // set rhoCp for fluid regions from thermofluid1
    Info<< "Adding to  rhoCp\n" << endl;
    rhoCpFluid.set
    (
      i,
      new volScalarField
      (
        IOobject
        (
          "rhoCp",
          runTime.timeName(),
          fluidRegions[i],
          IOobject::NO_READ,
          IOobject::NO_WRITE
        ),
        rhoFluid[i]*thermoFluid[i].Cp()
      )
    );
    // save for ddt term
    rhoCpFluid[i].oldTime();

    fvOptionsFluid.set
    (
        i,
        new fv::options
        (
            fluidRegions[i]
        )
    );
 
    Info<< "Creating field kinetic energy K\n" << endl;
    KFluid.set
    (
      i,
      new volScalarField
      (
        "K",
        0.5*magSqr(UFluid[i])
      )
    );
    
        Info<< "    Correcting initialMassFluid\n" << endl;
    initialMassFluid[i] =
        fvc::domainIntegrate(rhoFluid[i]).value();

    const dictionary& pimpleDict = 
        fluidRegions[i].solutionDict().subDict("PIMPLE");
    
    pimpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);
    // set correctPhi boolean here itself
    correctPhiFluid.set(i,pimpleDict.getOrDefault("correctPhi",fluidRegions[i].dynamic()));
        


    // create AlphaFluxes
    IOobject alphaPhi10Header
    (
        "alphaPhi10",
        runTime.timeName(),
        fluidRegions[i],
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    const bool alphaRestart =
    alphaPhi10Header.typeHeaderOk<surfaceScalarField>(true);

    // MULES flux from previous time-step
    alphaPhiFluid.set
    (
        i,
        new surfaceScalarField
        (
            alphaPhi10Header,
            phiFluid[i]*fvc::interpolate(thermoFluid[i].alpha1())
        )
    );
    
    alphaPhiUnFluid.set
    (
      i,
      new surfaceScalarField
      (
        IOobject
        (
          "alphaPhiUn",
          runTime.timeName(),
          fluidRegions[i],
          IOobject::NO_READ,
          IOobject::NO_WRITE
        ),
        fluidRegions[i],
        dimensionedScalar(phiFluid[i].dimensions(),Zero)
      )
    );
    kappaEffFluid.set
    (
      i,
      new volScalarField
      (
        IOobject
        (
          "kappaEff",
          runTime.timeName(),
          fluidRegions[i],
          IOobject::NO_READ,
          IOobject::NO_WRITE
        ),
        thermoFluid[i].kappa()
      )
    );


  }
