/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

Description
    Calculates and outputs the mean and maximum Courant Numbers for the fluid
    regions

\*---------------------------------------------------------------------------*/

#ifndef multiCourantNo_H
#define multiCourantNo_H

#include "fvMesh.H"
#include "multiphaseSystem.H"
// #include "mutiphaseSystems.H"

// add the function directly to namespace Foam so setMultiRegionDeltaT has access to it
namespace Foam
{
    //- returns alphaCoNum and effective fluid vel courant number in a pair 
    //- <CoNum, alphaCoNum> 
    Pair<scalar> multiCourantNoPair
    (
        const fvMesh& mesh,
        const Time& runTime,
        const surfaceScalarField& phi,
        const multiphaseSystem& fluid
    );
    
    // could pass just the phase system arg and init phi() from it in the function

    Field<scalar> multiCourantNo
    (
        const fvMesh& mesh,
        const Time& runTime,
        const surfaceScalarField& phi,
        const multiphaseSystem& fluid
    );
}

#endif

// ************************************************************************* //
