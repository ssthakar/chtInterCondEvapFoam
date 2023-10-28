    dynamicFvMesh& mesh = fluidRegions[i];
    surfaceScalarField& phi = phiFluid[i];// changed
    volScalarField& p_rgh = p_rghFluid[i];
    twoPhaseMixtureEThermo& thermo = thermoFluid[i]; // added 
    temperaturePhaseChangeTwoPhaseMixture&  mixture = mixtureFluid[i]; // added
    interfaceProperties& interface = interfaceFluid[i]; //added
    dimensionedScalar& rho1 = rho1Fluid[i]; //added
    dimensionedScalar& rho2 = rho2Fluid[i]; //added
    dimensionedScalar& Prt = PrtFluid[i]; // added
    volScalarField& alpha1 = thermoFluid[i].alpha1();
    volScalarField& alpha2 = thermoFluid[i].alpha2();
    surfaceScalarField& alphaPhi10 = alphaPhiFluid[i]; // added
    surfaceScalarField& alphaPhiUn = alphaPhiUnFluid[i]; // added
    tmp<surfaceScalarField> talphaPhi1Corr0; //added
    volScalarField& rho = rhoFluid[i];
    volScalarField& p = thermoFluid[i].p(); // changed
    volVectorField& U = UFluid[i];
    volScalarField& T = thermoFluid[i].T(); // should i put it in createFields??
    volScalarField& rhoCp = rhoCpFluid[i];
    surfaceScalarField& rhoPhi = rhoPhiFluid[i];
    volScalarField& gh = ghFluid[i];
    surfaceScalarField& ghf = ghfFluid[i];
    const dimensionedScalar& ghRef = ghRefFluid[i];
    volScalarField& K = KFluid[i];
    fv::options &fvOptions = fvOptionsFluid[i];
    volScalarField& kappaEff = kappaEffFluid[i];
    incompressible::turbulenceModel & turbulence = turbFluid[i];// changed
    scalar& cumulativeContErr = cumulativeContErrFluid[i];
    //bool& correctPhi = correctPhiFluid[i];
    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );
    rho.oldTime();
    bool frozenFlow = frozenFlowFluid[i];

    const label pRefCell = pRefCellFluid[i];
    const scalar pRefValue = pRefValueFluid[i];

    // PIMPLE control
    pimpleControl& pimple = pimpleFluid[i];

    bool correctPhi
    (
        pimple.dict().getOrDefault("correctPhi", mesh.dynamic())
    );
    
    bool checkMeshCourantNo
    (
        pimple.dict().getOrDefault("checkMeshCourantNo", false)
    );

    bool moveMeshOuterCorrectors
    (
        pimple.dict().getOrDefault("moveMeshOuterCorrectors", false)
    );
    // this needs to be changed
    volScalarField rAU
    (
      IOobject
      (
        "rAU",
        runTime.timeName(),
        fluidRegions[i],
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      fluidRegions[i],
      dimensionedScalar("rAU", dimTime/dimDensity, 1) 
    );
    
    IOMRFZoneList& MRF = MRFfluid[i]; // added
    // keep this the same as we dont' need to change anything here
    autoPtr<surfaceVectorField> Uf;

    if (fluidRegions[i].dynamic())
    {
        Info<< "Constructing face velocity Uf\n" << endl;

        Uf.reset
        (
            new surfaceVectorField
            (
                IOobject
                (
                    "Uf",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(UFluid[i])
            )
        );
    }


