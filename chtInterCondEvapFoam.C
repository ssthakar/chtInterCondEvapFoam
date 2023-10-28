/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    chtMultiRegionTwoPhaseEulerFoam

Group
    grpHeatTransferSolvers

Description
    Transient solver for buoyant, turbulent fluid flow and solid heat
    conduction with conjugate heat transfer between solid and fluid regions.

    It solves a two-phase Euler approach on the fluid region.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H" //add
#include "EulerDdtScheme.H" //add
#include "localEulerDdtScheme.H" //add 
#include "CrankNicolsonDdtScheme.H" //add
#include "dynamicFvMesh.H" // add 
#include "pimpleControl.H"
#include "interfaceProperties.H"
#include "regionProperties.H"
#include "multiCourantNo.H" // add 
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "radiationModel.H" // only solid implementation, fluid region does not use this 
#include "fvOptions.H"
#include "coordinateSystem.H"
#include "loopControl.H"
#include "CorrectPhi.H"
#include "subCycle.H" // add 
#include "temperaturePhaseChangeTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "turbulenceModel.H"
#include "twoPhaseMixtureEThermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for buoyant, turbulent two phase fluid flow and"
        "solid heat conduction with conjugate heat transfer "
        "between solid and fluid regions."
    );

    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    //debug switch for printing out info statements in code 
    #include "custom_debug.H"   
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "readSolidTimeControls.H"
    #include "multiPhaseMultiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"
    #include "setInitialMultiRegionDeltaT.H"

    while (runTime.run())
    {
        #include "readTimeControls.H" // maxCo declared here, reads from controlDict
        #include "readSolidTimeControls.H"
        #include "readPIMPLEControls.H"
        #include "multiPhaseMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (nOuterCorr != 1)
        {
            forAll(fluidRegions, i)
            {
                // #include "storeOldFluidFields.H"
            }
        }

        // --- PIMPLE loop same as while(pimple.loop())
        for (int oCorr=0; oCorr<nOuterCorr; ++oCorr)
        {   
            if(debug)
            {
              Info << "outer corrector val:  "<<oCorr<<"  "<<endl;
            } 
            const bool finalIter = (oCorr == nOuterCorr-1);
            forAll(fluidRegions, i)
            {
                Info<< "\nSolving for fluid region "
                    << fluidRegions[i].name() << endl;
                #include "setRegionFluidFields.H"
                #include "readFluidMultiRegionPIMPLEControls.H"
                #include "solveFluid.H"
            }
            forAll(solidRegions, i)
            {
                Info<< "\nSolving for solid region "
                    << solidRegions[i].name() << endl;
                #include "setRegionSolidFields.H"
                #include "readSolidMultiRegionPIMPLEControls.H"
                #include "solveSolid.H"
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
