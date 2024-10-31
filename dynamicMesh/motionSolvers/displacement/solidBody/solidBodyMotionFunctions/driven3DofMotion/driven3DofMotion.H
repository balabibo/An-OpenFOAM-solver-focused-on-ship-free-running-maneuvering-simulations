/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

Class
    Foam::solidBodyMotionFunctions::driven3DofMotion

Description
        
    Solidbodymotion which follows the motions of cOfG, including translation in X, Y axis and rotation in Z axis.


SourceFiles
    driven3DofMotion.C

\*---------------------------------------------------------------------------*/

#ifndef driven3DofMotion_H
#define driven3DofMotion_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"
#include "Function1.H"
#include "autoPtr.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class driven3DofMotion Declaration
\*---------------------------------------------------------------------------*/

class driven3DofMotion
:
    public solidBodyMotionFunction
{
    // Private data

        //- Name of the meshObject to dum CofG displacement
        word cOfGdisplacement_;

        //- Name of the meshObject to dum CofG rotation
        // word cOfGrotation_;

        //- Uniform vector to follow
        uniformDimensionedVectorField CofGdisp_;

        //- Uniform vector which stores cOfR ([0] [1]) and rotation angle ([2])
        uniformDimensionedVectorField CofGro_;     

    // Private Member Functions

        //- No copy construct
        driven3DofMotion(const driven3DofMotion&) = delete;

        //- No copy assignment
        void operator=(const driven3DofMotion&) = delete;


public:

    //- Runtime type information
    TypeName("driven3DofMotion");


    // Constructors

        //- Construct from components
        driven3DofMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new driven3DofMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~driven3DofMotion() = default;


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
