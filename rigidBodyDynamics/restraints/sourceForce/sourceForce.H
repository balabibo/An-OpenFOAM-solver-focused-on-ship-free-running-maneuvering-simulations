/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenFOAM Foundation
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
    Foam::RBD::restraints::sourceForce

Group
    grpRigidBodyDynamicsRestraints

Description
    Enforce the oumSource effects on prescribed rigid body

Usage
    Example applying a oumSource force to the hull:
    \verbatim
    restraints
    {
        sourceforce
        {
            type        sourceForce;
            body        hull;
            source      disk1;
            startTime   0;
        }
    }
    \endverbatim

SourceFiles
    sourceForce.C

\*---------------------------------------------------------------------------*/

#ifndef RBD_restraints_sourceForce_H
#define RBD_restraints_sourceForce_H

#include "rigidBodyRestraint.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace restraints
{

/*---------------------------------------------------------------------------*\
                          Class sourceForce Declaration
\*---------------------------------------------------------------------------*/

class sourceForce
:
    public restraint
{
    // Private data

        //- Name of momentum source term
        word source_;

        //- Moments and forces provided by the momentum source
        mutable spatialVector sourceForce_;

        //- Restraint startTime 
        scalar startTime_;



public:

    //- Runtime type information
    TypeName("sourceForce");


    // Constructors

        //- Construct from components
        sourceForce
        (
            const word& name,
            const dictionary& dict,
            const rigidBodyModel& model
        );

        //- Construct and return a clone
        virtual autoPtr<restraint> clone() const
        {
            return autoPtr<restraint>
            (
                new sourceForce(*this)
            );
        }


    //- Destructor
    virtual ~sourceForce();


    // Member Functions

        //- Accumulate the retraint internal joint forces into the tau field and
        //  external forces into the fx field
        virtual void restrain
        (
            scalarField& tau,
            Field<spatialVector>& fx,
            const rigidBodyModelState& state
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& dict);

        //- Write
        virtual void write(Ostream&) const;

        //- update sourceForce
        void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace restraints
} // End namespace RBD
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
