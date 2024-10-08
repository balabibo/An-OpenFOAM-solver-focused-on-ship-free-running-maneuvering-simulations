/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2018-2020 OpenCFD Ltd
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
// #include "fvMesh.H"
// #include "fvMatrix.H"
// #include "volFields.H"
// #include "IOobject.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::HOSource::calc
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<vector>& eqn,
    const labelList& cells
)
{
  const vectorField& U = eqn.psi();
  scalar refRho = 0;
  scalar Va = 0;
  forAll(wakeCells_, i)
  {
    refRho = refRho + rho[wakeCells_[i]];
    Va += mag((U[wakeCells_[i]] & diskDir_)*diskDir_ - (vcofR_ & diskDir_)*diskDir_); //axial velocity
  }
  reduce(refRho, sumOp<scalar>());
  reduce(Va, sumOp<scalar>());
  refRho /=  returnReduce(wakeCells_.size(), sumOp<label>());
  Va /=  returnReduce(wakeCells_.size(), sumOp<label>());

  vectorField& Usource = eqn.source();
  const vectorField& meshPosi = mesh().C();
  const scalarField& meshVolu = mesh().V();

  scalar J = Va/(diskRPS_*2*diskR_ + VSMALL);

  scalar thrust = KT(J)*refRho*sqr(diskRPS_)*pow4(2*diskR_);
  scalar torque = KQ(J)*refRho*sqr(diskRPS_)*pow5(2*diskR_);
  scalar   Ax   = 105*thrust/8/(M_PI*diskThick_*(3*diskH_*diskR_ + 4*diskR_)*(diskR_ - diskH_*diskR_));
  scalar Atheta = 105*torque/8/(M_PI*diskThick_*diskR_*(3*diskH_*diskR_ + 4*diskR_)*(diskR_ - diskH_*diskR_));
  // Info<<nl<<"target thrust = "<<thrust<<nl
      // <<nl<<"actual torque = "<<torque<<nl
    // <<nl<<"Ax = "<<Ax<<nl
    // <<nl<<"Atheta = "<<Atheta<<endl;
  fx_ = spatialVector(Zero, Zero);
  Info<<nl<<"Velocity of Rigid body = "<<vcofR_<<endl;

  scalar Tdisk = 0;
  scalar Qdisk = 0;

  forAll(cells, i)
  {

  vector OQ = meshPosi[cells[i]] - diskOri_;
  vector TQ = OQ - (OQ & diskDir_)*diskDir_ ;
  vector OT = -1*fabs(OQ & diskDir_) * diskDir_;
  // vector OT = OQ & diskDir_ * diskDir_;
  vector tanQ = rotaDir_*(OT ^ TQ)/(mag(OT ^ TQ)+VSMALL);
  // Pout<<nl<<"uAxis = "<<U[cells[i]]<<endl;
  scalar Rratio = mag(TQ)/diskR_;
  scalar Rstar  = (Rratio - diskH_)/(1 - diskH_);
  scalar   fbx  = Ax*Rstar*sqrt(fabs(1 - Rstar))*meshVolu[cells[i]];
  scalar fbtheta= Atheta*Rstar*sqrt(fabs(1 - Rstar))/(Rstar*(1 - diskH_) + diskH_)*meshVolu[cells[i]];
  Usource[cells[i]] += fbx*diskDir_ + fbtheta*tanQ;
  Tdisk += fbx;
  Qdisk += mag(TQ)*fbtheta;
  vector moment = meshPosi[cells[i]] ^ Usource[cells[i]];
  fx_ += spatialVector(moment, Usource[cells[i]]);

  }
  reduce(fx_, sumOp<spatialVector>());
  reduce(Tdisk, sumOp<scalar>());
  reduce(Qdisk, sumOp<scalar>());
  
  // Info<<nl<<"fx = "<<fx_<<endl;
  if
  (
    mesh().time().foundObject<uniformDimensionedSymmTensorField>
    (
        actBody_ + "sourceForce"
    )
  )
  {
    auto& fxx =
        mesh().time().lookupObjectRef<uniformDimensionedSymmTensorField>
        (
            actBody_ + "sourceForce"
        );
    for(label i =0; i < 6; i++)    
      fxx.value()[i] = fx_[i];
  }


  Info<<"J = "<<J<<nl
      <<"refRho = "<<refRho<<nl
      <<"target thrust = "<<thrust<<nl
      <<"actual thrust = "<<Tdisk<<nl
      <<"target KT = "<<thrust/(refRho*diskRPS_*diskRPS_*pow(2*diskR_,4))<<nl
      <<"actual KT = "<<Tdisk/(refRho*diskRPS_*diskRPS_*pow(2*diskR_,4))<<nl
      <<"target torque = "<<torque<<nl
      <<"actual torque = "<<Qdisk<<nl
      <<"target 10KQ = "<<10*torque/(refRho*diskRPS_*diskRPS_*pow(2*diskR_,5))<<nl
      <<"actual 10KQ = "<<10*Qdisk/(refRho*diskRPS_*diskRPS_*pow(2*diskR_,5))<<nl
      << "external force = "<<fx_<<endl;
  
  scalar Fy = 0;
  if(Pstream::master())
  {
    Ostream& os = file();
    writeCurrentTime(os);
    os << tab << diskRPS_<< tab <<Tdisk<< tab << Tdisk/(refRho*diskRPS_*diskRPS_*pow(2*diskR_,4))<< tab << Qdisk << tab << 10*Qdisk/(refRho*diskRPS_*diskRPS_*pow(2*diskR_,5)) << tab <<Fy <<endl;      
  }                               
    
/*
    if (mesh().time().outputTime())
    {
        volVectorField momentumSourceField
        (
            IOobject
            (
                "momentumSourceField",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            vector::zero
        );
        forAll(cells, i)
            momentumSourceField[cells[i]] = Usource[cells[i]];
        momentumSourceField.write();
    }
  */
    
}




// ************************************************************************* //
