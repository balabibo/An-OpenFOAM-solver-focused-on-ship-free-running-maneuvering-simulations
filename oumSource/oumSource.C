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

#include "oumSource.H"
#include "geometricOneField.H"
//#include "cellSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(oumSource, 0);
    addToRunTimeSelectionTable(option, oumSource, dictionary);
}
}




// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
/*
void Foam::fv::oumSource::writeFileHeader(Ostream& os)
{
    writeFile::writeHeader(os, "Actuation disk source");
    writeFile::writeCommented(os, "Time");
    writeFile::writeCommented(os, "Uref");
    writeFile::writeCommented(os, "Cp");
    writeFile::writeCommented(os, "Ct");

    if (forceMethod_ == forceMethodType::FROUDE)
    {
        writeFile::writeCommented(os, "a");
        writeFile::writeCommented(os, "T");
    }
    else if (forceMethod_ == forceMethodType::VARIABLE_SCALING)
    {
        writeFile::writeCommented(os, "Udisk");
        writeFile::writeCommented(os, "CpStar");
        writeFile::writeCommented(os, "CtStar");
        writeFile::writeCommented(os, "T");
        writeFile::writeCommented(os, "P");
    }

    os  << endl;
}
*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::oumSource::oumSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(name, modelType, dict, mesh),
    writeFile(mesh, name, modelType, coeffs_),
    rotaDir_(coeffs_.getOrDefault<label>("rotationDir", 1)),
    chordCoef_(coeffs_.getOrDefault<scalarList>("chordCoeff", scalarList(6,0))),
    pitchCoef_(coeffs_.getOrDefault<scalarList>("pitchCoeff", scalarList(6,0))),
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
    bladeN_(coeffs_.getOrDefault<label>("bladeNumber", 4)),
    diskRPS_(coeffs_.getOrDefault<scalar>("diskRPS", 5.0)),
    refBody_(coeffs_.getOrDefault<word>("refBody", "none")),
    actBody_(coeffs_.getOrDefault<word>("actBody", "none")),
    cofR_(Zero), 
    diskV_(0.0),
    fx_(Zero, Zero)
{
    Info << tab << "- creating actuation disk zone: " << this->name() << endl;

    
    // select cells forming part of the bodyforce disk
    updateCells();

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
    // Info<<nl<<""<<ffx.value()<<endl;
    // // Create a dimensioned origin position to use with default OF mesh types
    // dimensionedVector x0 ("x0", dimLength, diskOri_);

    // // compute distance and normal vectors from origin of the disk to use for selection
    // scalarField R (mag(mesh.C() - x0));
    // vectorField rHat ((mesh.C() - x0) / mag(mesh.C() - x0));

    // // go over each cell in the grid and comapre it against selection criteria
    // for (label cellI = 0; cellI < mesh.C().size(); cellI++)
    // {
    //     // determine distance from cell centre to disk axis and along the normal direction
    //     scalar dNormal = R[cellI] * (rHat[cellI] & diskDir_);
    //     scalar dRadial = sqrt(pow(R[cellI], 2.) - pow(dNormal, 2.));

    //     // see if the cell is within tolerance from the centre of the disk
    //     if ((mag(dNormal) < diskThick_/2.) && (dRadial < diskR_) && (dRadial > diskH_*diskR_))
    //         cells_.append(cellI);
    // }

    // ===
    // Set volume information
    // diskV_ = 0.0;
    // forAll(cells_, i)
    //     diskV_ += mesh.V()[cells_[i]];
    // reduce(diskV_, sumOp<scalar>());
    
    fieldNames_.resize(1, "U");

    fv::option::resetApplied();
       
    Info << tab
         << "- selected " << returnReduce(cells_.size(), sumOp<label>())
         << " cell(s) with volume " << diskV_ << endl;
    
    Info << nl
         << "refBody "<<refBody_ <<nl
         << "actBody " << actBody_ <<nl
         << "chordCoef_ " << chordCoef_ << nl
         << "pitchCoef_ " << pitchCoef_ << nl
         << "diskR_ " << diskR_ << nl
         << "diskH_ " << diskH_ << nl
         << "diskDir_" << diskDir_ << nl
         << "diskOri_" << diskOri_ << nl
         << "diskThick_ "<< diskThick_ << nl
         << "bladeN_ " <<bladeN_ << nl
         << "diskRPS_ " << diskRPS_<<endl;
         
    // tensor  roaTen(coeffs_.getOrDefault<tensor>("roaTensor",tensor::I));
    // vector  roaCenter(coeffs_.getOrDefault<vector>("roaCenter",Zero));
    // // // vector  omegaR(coeffs_.getOrDefault<vector>("omega",Zero));
    // vector  tranS(coeffs_.getOrDefault<vector>("v",Zero));
    // spatialTransform x00(I, vector(0.178, 0, 0.3323));
    // spatialTransform body1(roaTen, roaCenter);
    // spatialTransform X(body1.inv() & x00);
    // vector tranP = X.transformPoint(tranS);
    // vector body2 = spatialTransform(body1.E().T(), roaCenter).transformPoint(tranS);
    
    // // // spatialVector velo(omegaR, tranS);
    
    // vector vAfter = body1.transformPoint(tranS);
    // // spatialVector veloAfter = body1 & velo;
    // // spatialVector veloTranAfter = body1.T() & velo;
    // // spatialVector veloTranA = body1.inv() & velo;
    // Info<<nl<<"**********************************"<<nl
    //     <<"diskOri = "<<tranS<<nl
    //     <<"transformation = "<<body1<<nl
    //     <<"transformation & diskOri = "<<vAfter<<endl
    //     <<"transformation.T() & diskOri  = "<<body2<<endl
    //     <<"body1.inv() & x00 = "<<X<<nl
    //     <<"transformation  diskOri "<<tranP<<endl;
    //     <<"transformation.inv() = "<<body1.inv()<<nl
    //     <<"transformation.inv() & velo = "<<veloTranA<<nl
    //     <<endl;
     
    //  tensor aa(vector::uniform(1), Zero, Zero);
    //  aa = tensor::I & aa;
    //  Info<<nl<<"aa = "<<aa<<nl
    //  <<"aa & I = " <<aa<<nl
    //  <<"I & aa = " <<aa<<endl;      

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::oumSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    update();
    if (diskV() > VSMALL)
    {
        Info<<"!!!!!!!!!运行了!!!!!!!!!!!!!!!!!!——1"<<endl;
        calc(geometricOneField(), geometricOneField(), eqn, cells_);
    }
    else
    {
      Info<<"!!!!!!!!!没运行!!!!!!!!!!!!!!!!!!"<<endl;
    }
}


void Foam::fv::oumSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    update();
    if (diskV() > VSMALL)
    {
        Info<<"!!!!!!!!!运行了!!!!!!!!!!!!!!!!!!——2"<<endl;
        calc(geometricOneField(), rho, eqn, cells_);
    }
        else
    {
      Info<<"!!!!!!!!!没运行!!!!!!!!!!!!!!!!!!"<<endl;
    }
}


void Foam::fv::oumSource::addSup
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
        Info<<"!!!!!!!!!运行了!!!!!!!!!!!!!!!!!!——3"<<endl;
        calc(alpha, rho, eqn, cells_);
    }
        else
    {
      Info<<"!!!!!!!!!没运行!!!!!!!!!!!!!!!!!!"<<endl;
    }
}


bool Foam::fv::oumSource::read(const dictionary& dict)
{
    if (fv::option::read(dict))
    {
        dict.readIfPresent("rotationDir", rotaDir_);
        dict.readIfPresent("diskDir", diskDir_);
        diskDir_.normalise();
        if (mag(diskDir_) < VSMALL)
        {
            FatalIOErrorInFunction(dict)
                << "Actuator disk surface-normal vector is zero: "
                << "diskDir = " << diskDir_
                << exit(FatalIOError);
        }

        return true;
    }

    return false;
}

void Foam::fv::oumSource::update()
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
            refBody_ + "RotationCenter"
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
        const auto& constcCofR =
        mesh().lookupObject<uniformDimensionedVectorField>
        (
            refBody_ + "RotationCenter"
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

        cofR_ = constcCofR.value();
               
        spatialTransform X(consttCofR.value(), constvCofR.value());

        vector p = X.transformPoint(p_);
        diskOri_ = X.transformPoint(oriDiskOri_);
        diskDir_ = p -diskOri_;
        Info<<nl<<"momentum source position has been updated!"<<nl
            <<"transform = "<<X<<nl
            <<"p = "<<p<<nl
            <<"diskOri_ = "<<diskOri_<<nl
            <<"diskDir_ = "<<diskDir_<<endl;
    }
    cells_ = updateCells();    

    Info << tab
        << "- selected " << returnReduce(cells_.size(), sumOp<label>())
        << " cell(s) with volume " << diskV_ << endl;
    
    Info << nl
         << "refBody "<<refBody_ <<nl
         << "actBody " << actBody_ <<nl
         << "chordCoef_ " << chordCoef_ << nl
         << "pitchCoef_ " << pitchCoef_ << nl
         << "diskR_ " << diskR_ << nl
         << "diskH_ " << diskH_ << nl
         << "diskDir_" << diskDir_ << nl
         << "diskOri_" << diskOri_ << nl
         << "diskThick_ "<< diskThick_ << nl
         << "bladeN_ " <<bladeN_ << nl
         << "diskRPS_ " << diskRPS_<<endl;   
}

DynamicList<label> Foam::fv::oumSource::updateCells()
{
    DynamicList<label> cells;
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
        scalar dRadial = sqrt(pow(R[cellI], 2.) - pow(dNormal, 2.));

        // see if the cell is within tolerance from the centre of the disk
        if ((mag(dNormal) < diskThick_/2.) && (dRadial < diskR_) && (dRadial > diskH_*diskR_))
            cells.append(cellI);
    }

    // Set volume information
    diskV_ = 0.0;
    forAll(cells, i)
        diskV_ += mesh().V()[cells[i]];
    reduce(diskV_, sumOp<scalar>());

    Info<<nl<<"momentum source cells have been updated!"<<endl;
    return cells;
}


// ************************************************************************* //
