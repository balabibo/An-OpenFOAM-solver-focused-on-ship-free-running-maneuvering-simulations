/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "driven2DofMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "dimensionedVector.H"
//#include "rigidBodyMotion.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(driven2DofMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        driven2DofMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::driven2DofMotion::driven2DofMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    cOfGdisplacement_(SBMFCoeffs.get<word>("cOfGdisplacement")),
    CofGdisp_
    (
        IOobject
        (
            cOfGdisplacement_,
            time_.timeName(),
            "uniform",
            time_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dimensionedVector(dimless, Zero)
    ),
    CofGro_
    (
        IOobject
        (
            cOfGdisplacement_ + "3Dof",
            time_.timeName(),
            "uniform",
            time_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dimensionedVector(dimless, Zero)
    )
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::driven2DofMotion::transformation() const
{
    vector disp = CofGdisp_.value();
    disp[2]=0; // Only move in x and y direction

    DebugInFunction << "displacement  :" << CofGdisp_.value() << endl;
    quaternion R(1);
    septernion TR(septernion(-disp)*R); 

    DebugInFunction << "Time = " << time_.value()
                    << " transformation: " << TR << endl;
   
   vector axis(0, 0, 1); // The rotation axis Z
   vector origin = CofGro_.value();
   origin[2] = 0; // The center of rotation where 0 is in z direction
   quaternion RR(axis, 0.0);
   septernion TRR(septernion(-origin)*RR*septernion(origin));

   TR = TRR*TR; //The order of transformation should be from right to left, herein, the translation is conducted firstly, and the rotation is following.
    return TR;
}


bool Foam::solidBodyMotionFunctions::driven2DofMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    return true;
}


// ************************************************************************* //
