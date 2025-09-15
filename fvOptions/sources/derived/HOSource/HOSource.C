/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2018-2022 OpenCFD Ltd
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

#include "HOSource.H"
#include "geometricOneField.H"
//#include "cellSet.H"
#include "addToRunTimeSelectionTable.H"
#include "cellCellStencilObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(HOSource, 0);
    addToRunTimeSelectionTable(option, HOSource, dictionary);
}
}




// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::HOSource::HOSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(name, modelType, dict, mesh),
    writeFile(mesh, name, modelType, coeffs_),
    rotaDir_(1),
    KTCoef_(coeffs_.getOrDefault<scalarList>("KTCoeff", scalarList(6,0))),
    KQCoef_(coeffs_.getOrDefault<scalarList>("KQCoeff", scalarList(6,0))),
    diskR_
    (
        coeffs_.getCheck<scalar>
        (
            "diskRadius",
            scalarMinMax::ge(VSMALL)
        )
    ),
    diskH_
    (
        coeffs_.getCheck<scalar>
        (
            "diskHub",
            scalarMinMax::ge(VSMALL)
        )
    ),
    diskDir_
    (
        coeffs_.getCheck<vector>
        (
            "diskDir",
            [&](const vector& vec){ return mag(vec) > VSMALL; }
        ).normalise()
    ),
    oriDiskOri_(coeffs_.getOrDefault<vector>("diskOrigin", vector::uniform(0))),
    diskOri_(oriDiskOri_),
    p_(diskDir_ + oriDiskOri_),
    diskThick_
    (
        coeffs_.getCheck<scalar>
        (
            "diskThickness",
            scalarMinMax::ge(VSMALL)
        )
    ),
    diskRPS_(coeffs_.getOrDefault<scalar>("diskRPS", 5.0)),
    refBody_(coeffs_.getOrDefault<word>("refBody", "none")),
    actBody_(coeffs_.getOrDefault<word>("actBody", "none")),
    vcofR_(Zero), 
    diskV_(0.0),
    fx_(Zero, Zero)
{
    if(coeffs_.getOrDefault<word>("rotationDir", "right") == "left")
    {
        rotaDir_ = -1;
    }

    Info<< tab << "- creating actuation disk zone: " << this->name() 
        <<nl << tab <<  "- rotation direction: " << rotaDir_
        << endl;

    // select cells forming part of the bodyforce disk
    update();

    static uniformDimensionedSymmTensorField ffx
    (
        IOobject
        (
            actBody_ + "sourceForce",
            mesh.time().timeName(),
            "uniform",
            mesh.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dimensionedSymmTensor(dimless, Zero)
    );
    
    fieldNames_.resize(1, "U");

    fv::option::resetApplied();
       
    Info << tab
         << "- selected " << returnReduce(wakeCells_.size(), sumOp<label>())
         << " wake field cell(s) with volume " << diskV_ << endl
         << "- selected " << returnReduce(cells_.size(), sumOp<label>())
         << " cell(s) with volume " << diskV_ << endl;
    
    Info << nl
         << "refBody "<<refBody_ <<nl
         << "actBody " << actBody_ <<nl
         << "diskR_ " << diskR_ << nl
         << "diskH_ " << diskH_ << nl
         << "diskDir_" << diskDir_ << nl
         << "diskOri_" << diskOri_ << nl
         << "diskThick_ "<< diskThick_ << nl
         << "diskRPS_ " << diskRPS_<<endl;
        
    Ostream& os = file();
        os<<"Time"<<tab<<"diskRPS"<<tab<<"thrust"<<tab << "KT" <<tab<<"torque"<<tab<<"10KQ"<<tab<<"Fy"<<tab<<endl;

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::HOSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    update();
    if (diskV() > VSMALL)
    {
        Info<<"!!!!!!!!!Running!!!!!!!!!!!!!!!!!!——1"<<endl;
        calc(geometricOneField(), geometricOneField(), eqn, cells_);
    }
    else
    {
      Info<<"!!!!!!!!!diskV < VSMALL-Stop-!!!!!!!!!!!!!!!!!!"<<endl;
    }
}


void Foam::fv::HOSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    update();
    if (diskV() > VSMALL)
    {
        Info<<"!!!!!!!!!Running!!!!!!!!!!!!!!!!!!——2"<<endl;
        calc(geometricOneField(), rho, eqn, cells_);
    }
        else
    {
      Info<<"!!!!!!!!!diskV < VSMALL-Stop-!!!!!!!!!!!!!!!!!!"<<endl;
    }
}


void Foam::fv::HOSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    update();
    if (diskV() > VSMALL)
    {
        Info<<"!!!!!!!!!Running!!!!!!!!!!!!!!!!!!——3"<<endl;
        calc(alpha, rho, eqn, cells_);
    }
        else
    {
      Info<<"!!!!!!!!!diskV < VSMALL-Stop-!!!!!!!!!!!!!!!!!!"<<endl;
    }
}


bool Foam::fv::HOSource::read(const dictionary& dict)
{
    if (fv::option::read(dict))
    {
        return true;
    }

    return false;
}

void Foam::fv::HOSource::update()
{
    if
    (
        mesh().foundObject<uniformDimensionedScalarField>
        (
            actBody_
        )
    )
    {
        const auto& constoutputV =
        mesh().lookupObject<uniformDimensionedScalarField>
        (
            actBody_
        );  
        diskRPS_ = constoutputV.value()/2/M_PI;
        Info<<nl<<"revolution per second has been updated!"<<endl;
    }

    if
    (
        mesh().foundObject<uniformDimensionedVectorField>
        (
            refBody_ + "velocity"
        )
        &&        
        mesh().foundObject<uniformDimensionedVectorField>
        (
            refBody_ + "RotationVector"
        )
        &&
        mesh().foundObject<uniformDimensionedTensorField>
        (
            refBody_ + "RotationTensor"
        )        
    )
    {
        const auto& constvCCofR =
        mesh().lookupObject<uniformDimensionedVectorField>
        (
            refBody_ + "velocity"
        );

        const auto& constvCofR =
        mesh().lookupObject<uniformDimensionedVectorField>
        (
            refBody_ + "RotationVector"
        );

        const auto& consttCofR =
        mesh().lookupObject<uniformDimensionedTensorField>
        (
            refBody_ + "RotationTensor"
        );

        vcofR_ = constvCCofR.value();
               
        spatialTransform X(consttCofR.value(), constvCofR.value());

        vector p = X.transformPoint(p_);
        diskOri_ = X.transformPoint(oriDiskOri_);
        diskDir_ = p -diskOri_;
        Info<<nl<<"momentum source position has been updated!"<<endl;
    }

    updateCells();
    // cells_ = updateCells();
    // wakeCells_ = updateWakeCells();    

    Info << nl
        << "- selected " << returnReduce(wakeCells_.size(), sumOp<label>())
        << " wake field cell(s) with volume " << wakeCellV_ << endl
        << "- selected " << returnReduce(cells_.size(), sumOp<label>())
        << " cell(s) with volume " << diskV_ << endl;
    
    Info << nl
         << "diskR_ " << diskR_ << nl
         << "diskDir_" << diskDir_ << nl
         << "diskOri_" << diskOri_ << nl
         << "diskRPS_ " << diskRPS_<<endl;   
}

void Foam::fv::HOSource::updateCells()
{
    DynamicList<label> cells;
    DynamicList<label> wakeCells;
    if(mesh().foundObject<labelIOList>("zoneID"))
    {
        const Foam::cellCellStencilObject& overlap = Foam::Stencil::New(mesh());
        const labelList& cellTypes = overlap.cellTypes();
        dimensionedVector x0 ("x0", dimLength, diskOri_);
        scalarField R (mag(mesh().C() - x0));
        vectorField rHat ((mesh().C() - x0) / (mag(mesh().C() - x0)));
        for (label cellI = 0; cellI < mesh().C().size(); cellI++)
        {
            if(cellTypes[cellI] < 0.5)
            {
                
                scalar dNormal = R[cellI] * (rHat[cellI] & diskDir_);
                scalar dRadial = sqrt(fabs(pow(R[cellI], 2.) - pow(dNormal, 2.))); 
                if
                (
                    (fabs(dNormal) < diskThick_/2.) 
                    && (dRadial < diskR_) 
                    && (dRadial > diskH_*diskR_)
                )
                {   
                    cells.append(cellI);
                }

                else if
                (
                    (dNormal < 3*diskThick_/2.0)
                    && 
                    (dNormal > diskThick_) 
                    && 
                    (dRadial < diskR_) 
                )
                {
                    wakeCells.append(cellI);
                }
            }
        }
    }
    else
    {
        // Create a dimensioned origin position to use with default OF mesh types
        dimensionedVector x0 ("x0", dimLength, diskOri_);

        // compute distance and normal vectors from origin of the disk to use for selection
        scalarField R (mag(mesh().C() - x0));
        vectorField rHat ((mesh().C() - x0) / mag(mesh().C() - x0));

        // go over each cell in the grid and comapre it against selection criteria
        for (label cellI = 0; cellI < mesh().C().size(); cellI++)
        {
            // determine distance from cell centre to disk axis and along the normal direction
            scalar dNormal = R[cellI] * (rHat[cellI] & diskDir_);
            scalar dRadial = sqrt(fabs(pow(R[cellI], 2.) - pow(dNormal, 2.))); 

            // see if the cell is within tolerance from the centre of the disk
            if
            (
                (fabs(dNormal) < diskThick_/2.) 
                && (dRadial < diskR_) 
                && (dRadial > diskH_*diskR_)
            )
            {   
                cells.append(cellI);
            }
            else if
            (
                (dNormal < 4*diskThick_/2.0)
                && 
                (dNormal > 3*diskThick_/2.0) 
                && 
                (dRadial < diskR_) 
            )
            {
                wakeCells.append(cellI);
            }
        }
    }
    // Set volume information
    diskV_ = 0.0;
    forAll(cells, i)
        diskV_ += mesh().V()[cells[i]];
    reduce(diskV_, sumOp<scalar>());

    wakeCellV_ = 0.0;
    forAll(wakeCells, i)
        wakeCellV_ += mesh().V()[wakeCells[i]];
    reduce(wakeCellV_, sumOp<scalar>());

    cells_ = cells;
    wakeCells_ = wakeCells;
    Info<<nl<<"momentum source cells have been updated!"<<endl;

}


// ************************************************************************* //
