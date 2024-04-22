/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::rigidBodyMeshMotion

Group
    grpMeshMotionSolvers

Description
    Rigid-body mesh motion solver for fvMesh.

    Applies septernion interpolation of movement as function of distance to the
    object surface.

SourceFiles
    rigidBodyMeshMotion.C

\*---------------------------------------------------------------------------*/

#ifndef rigidBodyMeshMotion_H
#define rigidBodyMeshMotion_H

#include "displacementMotionSolver.H"
#include "rigidBodyMotion.H"
#include "ramp.H"
#include "maneuveringOutput.H"
#include "writeFile.H"
#include "forceMultiphase.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                     Class rigidBodyMeshMotion Declaration
\*---------------------------------------------------------------------------*/

class rigidBodyMeshMotion
:
    public displacementMotionSolver
{
    //- Class containing the patches and point motion weighting for each body
    class bodyMesh
    {
        //- Name of the body
        const word name_;

        //- ID of the body in the RBD::rigidBodyMotion
        const label bodyID_;

        //- List of mesh patches associated with this body
        const wordRes patches_;

        //- Patches to integrate forces
        const labelHashSet patchSet_;

        //- Inner morphing distance (limit of solid-body region)
        const scalar di_;

        //- Outer morphing distance (limit of linear interpolation region)
        const scalar do_;

        //- Current interpolation weight
        //  (1 at patches and within di_, 0 at do_ and beyond)
        pointScalarField weight_;

        //- Indicate the parent type of a body
        //word parentType_;

    public:

        friend class rigidBodyMeshMotion;

        bodyMesh
        (
            const polyMesh& mesh,
            const word& name,
            const label bodyID,
            const dictionary& dict
        );
    };


    // Private data

        //- Rigid-body model
        RBD::rigidBodyMotion model_;

        //- List of the bodyMeshes containing the patches and point motion
        //  weighting for each body
        PtrList<bodyMesh> bodyMeshes_;

        //  List of moving controllers which are used to change motions of 
        //  several rigid bodies
        PtrList<maneuveringOutput> maneuveringOutput_;

        //  List of uniformDimensionedScalarField for moving controllers 
        //  which is used to transmit rps for body-force source term
        // static PtrList<uniformDimensionedScalarField> outputValue_;

        //- Test-mode in which only the gravitational body-force is applied
        bool test_;

        //- Reference density required by the forces object for
        //  incompressible calculations, required if rho == rhoInf
        scalar rhoInf_;

        //- Name of density field, optional unless used for an
        //  incompressible simulation, when this needs to be specified
        //  as rhoInf
        word rhoName_;

        //- Ramp the forces according to the specified function and period
        autoPtr<Function1<scalar>> ramp_;

        //-Ramp time for rigid body motion, avoiding the initial crush of simulations
        scalar rampTime_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- Name of the uniformVectorField for CofG displacement
        word cOfGdisplacement_;

        // //- Body Id for the body to write CofG displacement
        // label bodyIdCofG_;
        //- Name of the body to write CofG displacement
        word bodyIdCofG_;
        
        //-output rigid body motions for post-processing 
        PtrList<functionObjects::writeFile> bodyMotionsOs_;

        //-output forces and moments acted on rigid body for post-processing 
        PtrList<functionObjects::writeFile> bodyDynamicsOs_;

        // //-output forces acted on rigid body for post-processing 
        // PtrList<functionObjects::writeFile> bodyForceOs_;

        // //-output moments acted on rigid body for post-processing 
        // PtrList<functionObjects::writeFile> bodyMomentOs_;

        //-forces and moments output file for rigid body on its local reference frame
        PtrList<functionObjects::forceMultiphase> forceFile_;

    // Private Member Functions

        //- No copy construct
        rigidBodyMeshMotion(const rigidBodyMeshMotion&) = delete;

        //- No copy assignment
        void operator=(const rigidBodyMeshMotion&) = delete;


public:

    //- Runtime type information
    TypeName("rigidBodyMotion");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        rigidBodyMeshMotion
        (
            const polyMesh&,
            const IOdictionary& dict
        );


    //- Destructor
    ~rigidBodyMeshMotion() = default;


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- acquire the rigid body's sailing velocity along x axis and yaw angle along z axis 
        const vector2D acquireInput(const word& refBody);

        //- alter the motions of objective rigid body
        void actControl(const word& actBody, const word& refBody, const label& controlType, const scalar& output);

        //- Solve for motion
        virtual void solve();

        //- Write state using stream options
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc
        ) const;

        //- Read dynamicMeshDict dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
