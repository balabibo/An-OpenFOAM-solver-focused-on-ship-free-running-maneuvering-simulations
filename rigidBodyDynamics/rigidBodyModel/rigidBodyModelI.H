/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::Time& Foam::RBD::rigidBodyModel::time() const
{
    return time_;
}


inline Foam::label Foam::RBD::rigidBodyModel::nBodies() const
{
    return bodies_.size();
}


inline Foam::PtrList<Foam::RBD::rigidBody>
Foam::RBD::rigidBodyModel::bodies() const
{
    return bodies_;
}

inline const Foam::DynamicList<Foam::label>& 
Foam::RBD::rigidBodyModel::countJoint() const
{
    return countJoint_;
}

inline const Foam::DynamicList<Foam::label>&
Foam::RBD::rigidBodyModel::lambda() const
{
    return lambda_;
}


inline const Foam::PtrList<Foam::RBD::joint>&
Foam::RBD::rigidBodyModel::joints() const
{
    return joints_;
}


inline Foam::label Foam::RBD::rigidBodyModel::nDoF() const
{
    return nDoF_;
}

inline const Foam::spatialVector Foam::RBD::rigidBodyModel::detailDof
(
    const label i
) const
{
    return detailDof_[i];
}

inline bool Foam::RBD::rigidBodyModel::unitQuaternions() const
{
    return unitQuaternions_;
}


inline const Foam::vector& Foam::RBD::rigidBodyModel::g() const
{
    return g_;
}


inline Foam::vector& Foam::RBD::rigidBodyModel::g()
{
    return g_;
}


inline const Foam::word& Foam::RBD::rigidBodyModel::name
(
    const label bodyID
) const
{
    if (merged(bodyID))
    {
        return mergedBody(bodyID).name();
    }
    else
    {
        return bodies_[bodyID].name();
    }
}


inline const Foam::RBD::rigidBodyInertia&
Foam::RBD::rigidBodyModel::I(const label i) const
{
    return bodies_[i];
}


// inline const Foam::spatialVector&
// Foam::RBD::rigidBodyModel::d(const label i) const
// {
//     return d_[i];
// }

inline const Foam::spatialVector&
Foam::RBD::rigidBodyModel::v(const label i) const
{
    return v_[i];
}

inline const Foam::spatialVector&
Foam::RBD::rigidBodyModel::a(const label i) const
{
    return a_[i];
}

inline bool Foam::RBD::rigidBodyModel::merged(label bodyID) const
{
    return bodyID < 0;
}


inline Foam::label Foam::RBD::rigidBodyModel::master(label bodyID) const
{
    if (bodyID < 0)
    {
        return mergedBody(bodyID).masterID();
    }
    else
    {
        return bodyID;
    }
}


inline Foam::label
Foam::RBD::rigidBodyModel::mergedBodyID(const label mergedBodyIndex) const
{
    return -1 - mergedBodyIndex;
}


inline Foam::label
Foam::RBD::rigidBodyModel::mergedBodyIndex(const label mergedBodyID) const
{
    return -1 - mergedBodyID;
}


inline const Foam::RBD::subBody&
Foam::RBD::rigidBodyModel::mergedBody(label mergedBodyID) const
{
    if (!merged(mergedBodyID))
    {
        FatalErrorInFunction
            << "Body " << mergedBodyID << " has not been merged"
            << abort(FatalError);
    }

    return mergedBodies_[mergedBodyIndex(mergedBodyID)];
}


inline Foam::label Foam::RBD::rigidBodyModel::bodyID(const word& name) const
{
    return bodyIDs_[name];
}


inline Foam::HashTable<Foam::label> Foam::RBD::rigidBodyModel::hbsize() const
{
    return bodyIDs_;
}


inline Foam::vector Foam::RBD::rigidBodyModel::masterPoint
(
    const label bodyID,
    const vector& p
) const
{
    if (merged(bodyID))
    {
        return
        (
            mergedBody(bodyID).masterXT().inv()
         && spatialVector(Zero, p)
        ).l();
    }
    else
    {
        return p;
    }
}


inline Foam::spatialVector Foam::RBD::rigidBodyModel::v
(
    const label bodyID,
    const vector& p
) const
{
    return
    (
        spatialTransform
        (
            X0_[master(bodyID)].E().T(),
            masterPoint(bodyID, p)
        )
      & v_[master(bodyID)]
    );
}


// ************************************************************************* //
