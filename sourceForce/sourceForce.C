/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "sourceForce.H"
#include "rigidBodyModel.H"
#include "rigidBodyModelState.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace restraints
{
    defineTypeNameAndDebug(sourceForce, 0);

    addToRunTimeSelectionTable
    (
        restraint,
        sourceForce,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::restraints::sourceForce::sourceForce
(
    const word& name,
    const dictionary& dict,
    const rigidBodyModel& model
)
:
    restraint(name, dict, model),
    source_(dict.getOrDefault<word>("source", "none")),
    sourceForce_(Zero, Zero),
    startTime_(dict.getOrDefault<scalar>("startTime", 0.0))
{
    read(dict);
    update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RBD::restraints::sourceForce::~sourceForce()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::restraints::sourceForce::restrain
(
    scalarField& tau,
    Field<spatialVector>& fx,
    const rigidBodyModelState& state
) const
{
    // Accumulate the force for the restrained body
    if(model_.time().timeOutputValue() > startTime_)
    {
        if
        (
        model_.time().foundObject<uniformDimensionedSymmTensorField>
        (
            source_ + "sourceForce"
        )
        )
        {
        auto& fxx =
            model_.time().lookupObjectRef<uniformDimensionedSymmTensorField>
            (
                source_ + "sourceForce"
            );
        for(label i =0; i < 6; i++)    
            sourceForce_[i] = fxx.value()[i];
        }

        if (model_.debug)
        {
            Info<< " source " << source_
                << " sourceForce " << sourceForce_
                << endl;
        }
        fx[bodyIndex_] += sourceForce_;
        Info<<nl<<"sourceForce has been acted on, sourceForce_ = "<<sourceForce_<<endl;
    }

}


bool Foam::RBD::restraints::sourceForce::read
(
    const dictionary& dict
)
{
    restraint::read(dict);

    coeffs_.readEntry("source", source_);
    coeffs_.readEntry("startTime", startTime_);

    // sourceForce_ = Function1<vector>::New("force", coeffs_);

    return true;
}


void Foam::RBD::restraints::sourceForce::write
(
    Ostream& os
) const
{
    restraint::write(os);

    os.writeEntry("source", source_);
    os.writeEntry("startTime", startTime_);

    // sourceForce_.writeData(os);
}


void Foam::RBD::restraints::sourceForce::update()
{

    if
    (
      model_.time().foundObject<uniformDimensionedSymmTensorField>
      (
          source_ + "sourceForce"
      )
    )
    {
      auto& fxx =
          model_.time().lookupObjectRef<uniformDimensionedSymmTensorField>
          (
              source_ + "sourceForce"
          );
      for(label i =0; i < 6; i++)    
        sourceForce_[i] = fxx.value()[i];
    }
}


// ************************************************************************* //
