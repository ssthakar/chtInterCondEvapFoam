
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

Global
    alphaCourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

\*---------------------------------------------------------------------------*/
// get the user specified courant numbers from controlDict

scalar maxAlphaCo
(
    runTime.controlDict().get<scalar>("maxAlphaCo")
);

scalar maxAlphaDdt
(
    runTime.controlDict().getOrDefault("maxAlphaDdt", GREAT)
);

scalar maxDiFluid
(
    runTime.controlDict().getOrDefault<scalar>("maxDiFluid", GREAT)
);

scalar CoNum = -GREAT;
scalar alphaCoNum = -GREAT;
//scalar meanAlphaCoNum = -GREAT;

scalar ddtAlphaNum = -GREAT;
scalar DiNumFluid = -GREAT;

// loop through all fluid regions 
forAll(fluidRegions, regionI)
{
  Field<scalar> CourantNos = multiCourantNo
  (
    fluidRegions[regionI],
    runTime,
    phaseSystemFluid[regionI].phi(), // get phi from phase system
    phaseSystemFluid[regionI] //phase system itself
  );
  CoNum = max
  (
      CourantNos[0],
      CoNum
  );

  alphaCoNum = max
  (
      CourantNos[1],
      alphaCoNum
  );
  ddtAlphaNum = max
  (
    CourantNos[2],
    ddtAlphaNum
    
  );
  DiNumFluid =max
  (
    CourantNos[3],
    DiNumFluid
  );
}

// ************************************************************************* //
