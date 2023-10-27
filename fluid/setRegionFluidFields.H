    dynamicFvMesh& mesh = fluidRegions[i];
    multiphaseSystem& fluid  = phaseSystemFluid[i];
    surfaceScalarField& phi = fluid.phi();
    volScalarField& p_rgh = p_rghFluid[i];
    volScalarField& rho = rhoFluid[i];
    volScalarField& p = pFluid[i];
    volVectorField& U = UFluid[i];
    volScalarField& T = TFluid[i];
    volScalarField& rhoCp = rhoCpFluid[i];
    surfaceScalarField& rhoPhi = fluid.rhoPhi();
    volScalarField& gh = ghFluid[i];
    surfaceScalarField& ghf = ghfFluid[i];
    const dimensionedScalar& ghRef = ghRefFluid[i];
    volScalarField& K = KFluid[i];
    fv::options &fvOptions = fvOptionsFluid[i];
    CompressibleTurbulenceModel<multiphaseSystem> & turbulence = turbFluid[i];
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


