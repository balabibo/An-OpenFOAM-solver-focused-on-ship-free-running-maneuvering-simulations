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

#include "oumSource.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "volFields.H"
// #include "IOobject.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::oumSource::calc
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<vector>& eqn,
    const labelList& cells
)
{
  const vectorField& U = eqn.psi();
  vectorField& Usource = eqn.source();
  const vectorField& meshPosi = mesh().C();
  const scalarField& meshVolu = mesh().V();
  scalar thrust = 0;
  scalar torque = 0;
  const scalar k1 = 1.07 - 1.05*(chordLength(0.7)*0.7) + 0.375*pow(chordLength(0.7)*0.7,2);
  const scalar cD = 0.01;
  fx_ = spatialVector(Zero, Zero);
  forAll(cells, i)
  {
    vector OQ = meshPosi[cells[i]] - diskOri_;
    vector TQ = OQ - (OQ & diskDir_)*diskDir_ ;
    vector OT = abs(OQ & diskDir_) * diskDir_;
    vector tanQ = (OT ^ TQ)/(mag(OT ^ TQ)+VSMALL);
    vector uAxis = (U[cells[i]] & diskDir_)*diskDir_; //axial velocity
    vector uTang = (U[cells[i]] & tanQ)*tanQ; //tangential velocity
    scalar diskR = mag(TQ);
    scalar relativeR = diskR/diskR_;
    scalar uTotal = sqrt(magSqr(uAxis)+pow(2*M_PI*diskRPS_*diskR - mag(uTang),2)); // resultant velocity
    //scalar uTotal = sqrt(magSqr(uAxis) + magSqr((2*M_PI*diskRPS_*diskDir_ ^ TQ) - uTang));
    scalar hydroPitch = atan(mag(uAxis)/(2*M_PI*diskRPS_*diskR - mag(uTang))); // hydrodyanmic pitch angle
    scalar attackAngle = atan(pitchValue(relativeR)/(2*M_PI)) - hydroPitch;//-- attack angle
    scalar cL = 2*M_PI*k1*sin(attackAngle);
    scalar fAxial = alpha[cells[i]]*rho[cells[i]]*meshVolu[cells[i]] * bladeN_*0.5*uTotal*uTotal*chordLength(relativeR)*(cL*cos(hydroPitch) - cD*sin(hydroPitch))/diskThick_/(2*M_PI);
    scalar fTang = alpha[cells[i]]*rho[cells[i]]*meshVolu[cells[i]] * bladeN_*0.5*uTotal*uTotal*chordLength(relativeR)*(cL*sin(hydroPitch) + cD*cos(hydroPitch))/diskThick_/(2*M_PI);
    thrust += fAxial;
    torque += fTang*diskR;
    Usource[cells[i]] = fAxial*diskDir_ - fTang*tanQ;
    vector moment = meshPosi[cells[i]] ^ Usource[cells[i]];
    fx_ += spatialVector(moment, Usource[cells[i]]);
  }
    reduce(thrust, sumOp<scalar>());
    reduce(torque, sumOp<scalar>());
    reduce(fx_, sumOp<spatialVector>());

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

    
    scalar refRho = 0;
    scalar refAlpha = 0;
    scalar Fy = 0;
    scalar Fz = 0;
    // scalar zB = 0.1;
    /*
      Info << "此前源项值： "<< Usource[cells[0]] << endl; // 此处语句会导致程序出错，但表面上看并无问题，放在循环中读取cell[0]也无问题
    */
    forAll(cells, i)
    {

      //  Info <<"i = "<< i << "  此前源项值： "<< Usource[cells[i]] << endl;
    //  if(fabs(meshPosi[cells[i]][0])<zB)
    //  {
    //    zB = fabs(meshPosi[cells[i]][0]);
    //  }
     

      refRho = refRho + rho[cells[i]];
      refAlpha = refAlpha + alpha[cells[i]];
      Fy += Usource[cells[i]][1];
      Fz += Usource[cells[i]][2];
    }
    //  Info<<"x坐标最小为： "<<zB<<endl;
     reduce(refRho, sumOp<scalar>());
     reduce(refAlpha, sumOp<scalar>());
     reduce(Fy, sumOp<scalar>());
     reduce(Fz, sumOp<scalar>());
     refRho /=  returnReduce(cells_.size(), sumOp<label>());
     refAlpha /= returnReduce(cells_.size(), sumOp<label>());
     
     
    //  Info << "现在源项值： " << Usource[cells[0]] << endl;
    //  Info<<"J = "<<refU/(diskRPS_*2*diskR_)<<endl
    Info<<"refRho = "<<refRho<<endl
        <<"refAlpha = "<<refAlpha<<endl
        <<"thrust = "<<thrust<<endl
        <<"KT:  "<<thrust/(refRho*diskRPS_*diskRPS_*pow(2*diskR_,4))<<endl
        <<"10KQ:  "<<10*torque/(refRho*diskRPS_*diskRPS_*pow(2*diskR_,5))<<endl
        << "Fy = "<<Fy<<endl
        << "Fz = "<<Fz<<endl
        << "external force = "<<fx_<<endl;
      if(thrust)
      {
        Ostream& os = file();
        writeCurrentTime(os);
        os << tab <<thrust<< tab << thrust/(refRho*diskRPS_*diskRPS_*pow(2*diskR_,4))<< tab <<Fy <<tab <<Fz <<endl;      
      }                               
    
    /*
    if
    (
        mesh_.time().timeOutputValue() >= writeFileStart_
     && mesh_.time().timeOutputValue() <= writeFileEnd_
    )
    {
        Ostream& os = file();
        writeCurrentTime(os);

        os  << Uref << tab << Cp << tab << Ct << tab << a << tab << T
            << endl;
    }
    */  
    
}




// ************************************************************************* //
