/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "rigidBodyMeshMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "pointPatchDist.H"
#include "pointConstraints.H"
#include "uniformDimensionedFields.H"
#include "forces.H"
#include "mathematicalConstants.H"
#include "maneuveringOutput.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rigidBodyMeshMotion, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        rigidBodyMeshMotion,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rigidBodyMeshMotion::bodyMesh::bodyMesh
(
    const polyMesh& mesh,
    const word& name,
    const label bodyID,
    const dictionary& dict
)
:
    name_(name),
    bodyID_(bodyID),
    patches_(dict.get<wordRes>("patches")),
    patchSet_(mesh.boundaryMesh().patchSet(patches_)),
    di_(dict.get<scalar>("innerDistance")),
    do_(dict.get<scalar>("outerDistance")),
    weight_
    (
        IOobject
        (
            name_ + ".motionScale",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, Zero)
    )
{}


Foam::rigidBodyMeshMotion::rigidBodyMeshMotion
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    model_
    (
        mesh.time(),
        coeffDict(),
        IOobject
        (
            "rigidBodyMotionState",
            mesh.time().timeName(),
            "uniform",
            mesh
        ).typeHeaderOk<IOdictionary>(true)
      ? IOdictionary
        (
            IOobject
            (
                "rigidBodyMotionState",
                mesh.time().timeName(),
                "uniform",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        )
      : coeffDict()
    ),
    test_(coeffDict().getOrDefault("test", false)),
    rhoInf_(1.0),
    rhoName_(coeffDict().getOrDefault<word>("rho", "rho")),
    ramp_
    (
        Function1<scalar>::NewIfPresent("ramp", coeffDict(), word::null, &mesh)
    ),
    rampTime_(coeffDict().getOrDefault<scalar>("rampTime", 0)),
    curTimeIndex_(-1),
    cOfGdisplacement_
    (
        coeffDict().getOrDefault<word>("cOfGdisplacement", "none")
    ),
    bodyIdCofG_(coeffDict().getOrDefault<word>("bodyIdCofG", "none"))
{
    if (rhoName_ == "rhoInf")
    {
        readEntry("rhoInf", rhoInf_);
    }

    dictionary maneuversDict = coeffDict().subOrEmptyDict("maneuvers");

    if
    (
      IOobject
      (
         "maneuvers",
         mesh.time().timeName(),
         "uniform",
         mesh
      ).typeHeaderOk<IOdictionary>(true)
    )
    {
        
        IOdictionary maneuvers
        (
            IOobject
            (
                "maneuvers",
                mesh.time().timeName(),
                "uniform",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        );
        maneuversDict = maneuvers;

        Info<<nl<<"*********************"<<nl<<"reading maneuvers dictionary!"<<nl<<maneuversDict<<endl;
    }

    static PtrList<uniformDimensionedScalarField> outputValue;
    static PtrList<uniformDimensionedVectorField> vectorofR;
    static PtrList<uniformDimensionedVectorField> vcenterofR;
    static PtrList<uniformDimensionedTensorField> transformCofR;
    for (const entry& dEntry : maneuversDict)
    {
        if(dEntry.isDict())
        {
            const dictionary& clDict = dEntry.dict();
            maneuveringOutput_.append
            (
                new maneuveringOutput
                (
                    mesh,
                    clDict                       
                )
            );
            if(clDict.get<word>("type") == "sailing")
            {
                bool oumFlag = true;
                for(label i = 0; i < model_.nBodies(); i++)
                {
                    if(model_.name(i) == clDict.get<word>("actBody"))
                    {
                        oumFlag = false;
                    }
                }
                if(oumFlag)
                {
                    oumFile_.append
                    (
                        new functionObjects::writeFile
                        (
                            mesh,
                            clDict.get<word>("actBody") + "Dynamics",
                            clDict.get<word>("actBody") + "Dynamics",
                            dict
                        )
                    );
                    oumName_.append
                    (
                        new word(clDict.get<word>("actBody"))
                    );
                }
                
                outputValue.append
                (
                    new uniformDimensionedScalarField
                    (
                        IOobject
                        (
                            clDict.get<word>("actBody"),
                            mesh.time().timeName(),
                            "uniform",
                            mesh,
                            IOobject::READ_IF_PRESENT,
                            IOobject::AUTO_WRITE
                        ),
                        dimensionedScalar(dimless, 0)
                    )

                );

                vcenterofR.append
                (
                    new uniformDimensionedVectorField
                    (
                        IOobject
                        (
                            clDict.get<word>("refBody") + "velocity",
                            mesh.time().timeName(),
                            "uniform",
                            mesh,
                            IOobject::READ_IF_PRESENT,
                            IOobject::NO_WRITE
                        ),
                        dimensionedVector(dimless, Zero)
                    )

                );

                vectorofR.append
                (
                    new uniformDimensionedVectorField
                    (
                        IOobject
                        (
                            clDict.get<word>("refBody") + "RotationVector",
                            mesh.time().timeName(),
                            "uniform",
                            mesh,
                            IOobject::READ_IF_PRESENT,
                            IOobject::AUTO_WRITE
                        ),
                        dimensionedVector(dimless, Zero)
                    )

                );

                transformCofR.append
                (
                    new uniformDimensionedTensorField
                    (
                        IOobject
                        (
                            clDict.get<word>("refBody") + "RotationTensor",
                            mesh.time().timeName(),
                            "uniform",
                            mesh,
                            IOobject::READ_IF_PRESENT,
                            IOobject::AUTO_WRITE
                        ),
                        dimensionedTensor(dimless, I)
                    )

                );

            }
        }
    }
    
    forAll(oumFile_, oi)
    {
        Ostream& os = oumFile_[oi].file();
        os<<"Time"<<tab;
        for(label i =1; i<= 6; i++)
            os<<"total"<<i<<tab;
        os<<endl;
    }

    const dictionary& bodiesDict = coeffDict().subDict("bodies");

    for (const entry& dEntry : bodiesDict)
    {
        const keyType& bodyName = dEntry.keyword();
        const dictionary& bodyDict = dEntry.dict();

        if (bodyDict.found("patches"))
        {
            const label bodyID = model_.bodyID(bodyName);
            //const word parentType(bodyDict.get<word>("parent"));
            if (bodyID == -1)
            {
                FatalErrorInFunction
                    << "Body " << bodyName
                    << " has been merged with another body"
                       " and cannot be assigned a set of patches"
                    << exit(FatalError);
            }

            bodyMeshes_.append
            (
                new bodyMesh
                (
                    mesh,
                    bodyName,
                    bodyID,
                    bodyDict
                )
            );
        }
    }

    // Calculate scaling factor everywhere for each meshed body
    forAll(bodyMeshes_, bi)
    {
        const pointMesh& pMesh = pointMesh::New(mesh);

        pointPatchDist pDist(pMesh, bodyMeshes_[bi].patchSet_, points0());

        pointScalarField& scale = bodyMeshes_[bi].weight_;

        // Scaling: 1 up to di then linear down to 0 at do away from patches
        scale.primitiveFieldRef() =
            min
            (
                max
                (
                    (bodyMeshes_[bi].do_ - pDist.primitiveField())
                   /(bodyMeshes_[bi].do_ - bodyMeshes_[bi].di_),
                    scalar(0)
                ),
                scalar(1)
            );

        // Convert the scale function to a cosine
        scale.primitiveFieldRef() =
            min
            (
                max
                (
                    0.5
                  - 0.5
                   *cos(scale.primitiveField()
                   *Foam::constant::mathematical::pi),
                    scalar(0)
                ),
                scalar(1)
            );

        pointConstraints::New(pMesh).constrain(scale);
        //scale.write();
    }



    forAll(bodyMeshes_, bi)
    {
        bodyMotionsOs_.append
        (
            new functionObjects::writeFile
            (
                mesh,
                bodyMeshes_[bi].name_ + "Motion",
                bodyMeshes_[bi].name_ + "Motion",
                dict
            )
        );
        Ostream& os = bodyMotionsOs_[bi].file();
        os<<"Time"<<tab;
        for(label i =1; i<= 6; i++)
            os<<"x"<<i<<tab;
        for(label i =1; i<= 6; i++)
            os<<"v"<<i<<tab;
        for(label i =1; i<= 6; i++)
            os<<"a"<<i<<tab;
        os<<endl;
    }
    
    forAll(bodyMeshes_, bi)
    {
        dictionary forcesDict;
        forcesDict.add("type", functionObjects::forceMultiphase::typeName);
        forcesDict.add("patches", bodyMeshes_[bi].patches_);
        forcesDict.add("phase", "water");
        forcesDict.add("CofR", vector::zero);
        // forcesDict.add("writeControl", "timeStep");
        // forcesDict.add("timeInterval", 1);

        forceFile_.append
        (
            new functionObjects::forceMultiphase
            (
                bodyMeshes_[bi].name_,
                db(),
                forcesDict
            )

        );

        bodyDynamicsOs_.append
        (
            new functionObjects::writeFile
            (
                mesh,
                bodyMeshes_[bi].name_ + "Dynamics",
                bodyMeshes_[bi].name_ + "Dynamics",
                dict
            )
        );
        Ostream& os = bodyDynamicsOs_[bi].file();
        os<<"Time"<<tab<<"wettedArea"<<tab;
        for(label i =1; i<= 6; i++)
            os<<"total"<<i<<tab;
        for(label i =1; i<= 6; i++)
            os<<"pressure"<<i<<tab;
        for(label i =1; i<= 6; i++)
            os<<"viscous"<<i<<tab;
        os<<endl;
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::rigidBodyMeshMotion::curPoints() const
{
    tmp<pointField> newPoints(points0() + pointDisplacement_.primitiveField());

    if (moveAllCells())
    {
        return newPoints;
    }
    else
    {
        tmp<pointField> ttransformedPts(new pointField(mesh().points()));
        pointField& transformedPts = ttransformedPts.ref();

        UIndirectList<point>(transformedPts, pointIDs()) =
            pointField(newPoints.ref(), pointIDs());

        return ttransformedPts;
    }
}

const vector2D Foam::rigidBodyMeshMotion::acquireInput(const word& refBody, const label& controlType)
{
    vector2D sailingState(Zero);
    sailingState[0] = this->model_.v(this->model_.bodyID(refBody), Zero)[3]; // sailing velocity along x axis
    // for (label i=0; i< this->model_.nDoF(); i++)
    // {
    //     Info<<this->model_.v(this->model_.bodyID(refBody), Zero)[2]<<nl
    //         <<this->model_.state().qDot()[i]
    //         <<endl;
    //     if (this->model_.v(this->model_.bodyID(refBody), Zero)[2] == this->model_.state().qDot()[i])
    //     {
    //         sailingState[1] = this->model_.state().q()[i]; // yaw angle along z axis
    //     }
    // }
    label positionDof = 0;
    forAll(bodyMeshes_, bi)
    {
        if(bodyMeshes_[bi].name_ == refBody)
        {
            positionDof = bi+1;
            sailingState[1] = model_.state().q()[model_.detailDof(positionDof)[controlType]];
        }
    }
    return sailingState;  

}


void Foam::rigidBodyMeshMotion::actControl(const word& actBody, const word& refBody, const label& controlType, const scalar& output)
{
    label positionDof = 0;
    forAll(bodyMeshes_, bi)
    {
        if(bodyMeshes_[bi].name_ == actBody)
        {
            positionDof = bi+1;
            // Info<<nl<<""<<model_.state().qDot()[model_.detailDof(positionDof)[controlType]]<<endl;
            model_.state().qDot()[model_.detailDof(positionDof)[controlType]] = output;

            Info<<nl<<"*****************"<<"outputControl = "<<output<<"*************"<<endl;

            model_.state().qDdot()[model_.detailDof(positionDof)[controlType]] 
            = 
            (model_.state().qDot()[model_.detailDof(positionDof)[controlType]] - model_.state0().qDot()[model_.detailDof(positionDof)[controlType]])/mesh().time().deltaTValue();

            model_.state().q()[model_.detailDof(positionDof)[controlType]]
            = 
            (model_.state().qDot()[model_.detailDof(positionDof)[controlType]] + model_.state0().qDot()[model_.detailDof(positionDof)[controlType]]) * mesh().time().deltaTValue()/2.0 + model_.state0().q()[model_.detailDof(positionDof)[controlType]];
        }
        
    }

    model_.forwardDynamicsCorrection(model_.state());

}

void Foam::rigidBodyMeshMotion::solve()
{
    const Time& t = mesh().time();

    if (mesh().nPoints() != points0().size())
    {
        FatalErrorInFunction
            << "The number of points in the mesh seems to have changed." << endl
            << "In constant/polyMesh there are " << points0().size()
            << " points; in the current mesh there are " << mesh().nPoints()
            << " points." << exit(FatalError);
    }

    // Store the motion state at the beginning of the time-step
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        model_.newTime();
        curTimeIndex_ = this->db().time().timeIndex();
    }

    const scalar ramp = (ramp_ ? ramp_->value(t.value()) : 1.0);

    if (t.foundObject<uniformDimensionedVectorField>("g"))
    {
        model_.g() =
            ramp*t.lookupObject<uniformDimensionedVectorField>("g").value();
    }

    vector oldPos(vector::uniform(GREAT));
    scalar oldW(0);
    if (bodyIdCofG_ != "none")
    {
        oldPos =model_.cCofR(model_.bodyID(bodyIdCofG_));
        oldW = model_.wCofR(model_.bodyID(bodyIdCofG_))[2];
    }

    if (test_)
    {
        const label nIter(coeffDict().get<label>("nIter"));

        for (label i=0; i<nIter; i++)
        {
            model_.solve
            (
                t.value(),
                t.deltaTValue(),
                scalarField(model_.nDoF(), Zero),
                Field<spatialVector>(model_.nBodies(), Zero)
            );
        }
    }
    else 
    {
        if (t.value() >= rampTime_)
        {
            const label nIter(coeffDict().getOrDefault<label>("nIter", 1));
            for (label i=0; i<nIter; i++)
            {
                Field<spatialVector> fx(model_.nBodies(), Zero);
                if(db().foundObject<volScalarField> ("p_rgh"))
                {
                    forAll(bodyMeshes_, bi)
                    {
                        const label bodyID = bodyMeshes_[bi].bodyID_;
                        dictionary forcesDict;
                        forcesDict.add("type", functionObjects::forceMultiphase::typeName);
                        forcesDict.add("patches", bodyMeshes_[bi].patches_);
                        // forcesDict.add("p", "p"); // p is the default value in fact
                        forcesDict.add("phase", "water");
                        forcesDict.add("CofR", vector::zero);

                        functionObjects::forceMultiphase f("flow on rigid body", db(), forcesDict);
                        
                        f.calcForcesMoments();

                        fx[bodyID] = ramp*spatialVector(f.momentEff(), f.forceEff());
                    }
                    
                }
                else
                {
                    forAll(bodyMeshes_, bi)
                    {
                        const label bodyID = bodyMeshes_[bi].bodyID_;
                        dictionary forcesDict;
                        forcesDict.add("type", functionObjects::forces::typeName);
                        forcesDict.add("patches", bodyMeshes_[bi].patches_);
                        forcesDict.add("rhoInf", rhoInf_);
                        forcesDict.add("rho", rhoName_);
                        forcesDict.add("CofR", vector::zero);

                        functionObjects::forces f("forces", db(), forcesDict);
                        f.calcForcesMoments();

                        fx[bodyID] = ramp*spatialVector(f.momentEff(), f.forceEff());
                    }
                }
                
                model_.solve
                (
                    t.value(),
                    t.deltaTValue(),
                    scalarField(model_.nDoF(), Zero),
                    fx
                );
            }
        
        }
        else
        {
            Info<<nl<<"Current time = "<<t.value()<<" is lower than ramp time = "<<rampTime_<<", its motion is restrained"<<endl;
        }

        forAll(maneuveringOutput_, bi)
        {
            const word refBody = maneuveringOutput_[bi].mInput().refBody();
            const word actBody = maneuveringOutput_[bi].mInput().actBody();
            const label controlType = maneuveringOutput_[bi].mInput().controlType();
            const scalar output = maneuveringOutput_[bi].output(this->acquireInput(refBody, controlType));
            
            bool meshFlag = true;
            for(label i = 0; i < model_.nBodies(); i++)
            {
                if(model_.name(i) == actBody)
                {
                    meshFlag = false;
                    //Info<<nl<<"*******************************"<<nl<<"Current Ref speed: "<<this->acquireInput(refBody)
                    //    <<nl<<"*******************************"<<endl;
                    this->actControl(actBody, refBody, controlType, output);
                }
            }

            if(meshFlag)
            {
                Info<<nl<<"creating uniformField for momentum source"<<endl;
                auto& outputV =
                    mesh().lookupObjectRef<uniformDimensionedScalarField>
                    (
                        actBody
                    );

                outputV.value() = output;

                auto& vCCofR =
                    mesh().lookupObjectRef<uniformDimensionedVectorField>
                    (
                        refBody + "velocity"
                    );
                // vCCofR.value() = model_.X0(model_.bodyID(refBody)).r();
                vCCofR.value() = model_.vCofR(model_.bodyID(refBody));

                spatialTransform X(model_.X0(model_.bodyID(refBody)).inv() & model_.X00(model_.bodyID(refBody)));

                auto& vCofR =
                    mesh().lookupObjectRef<uniformDimensionedVectorField>
                    (
                        refBody + "RotationVector"
                    );
                vCofR.value() = X.r();

                auto& tCofR =
                    mesh().lookupObjectRef<uniformDimensionedTensorField>
                    (
                        refBody + "RotationTensor"
                    );
                tCofR.value() = X.E();  
                
            }
        }

        vector presentPos(vector::uniform(GREAT));
        scalar presentW(0);
        if (bodyIdCofG_ != "none")
        {
            presentPos =model_.cCofR(model_.bodyID(bodyIdCofG_));
            presentW = model_.wCofR(model_.bodyID(bodyIdCofG_))[2];
        }

        if (cOfGdisplacement_ != "none")
        {
            if (bodyIdCofG_ != "none")
            {
                if
                (
                    db().time().foundObject<uniformDimensionedVectorField>
                    (
                        cOfGdisplacement_
                    )
                )
                {
                    auto& disp =
                        db().time().lookupObjectRef<uniformDimensionedVectorField>
                        (
                            cOfGdisplacement_
                        );

                    // disp.value() += model_.cCofR(bodyIdCofG_) - oldPos;
                    disp.value() += presentPos - oldPos;
                    // Info<<nl<<"presentPos = "<<presentPos<<nl
                    //     <<"oldPos = "<<oldPos<<nl
                    //     <<"disp.value() = "<<disp.value()<<endl;
                }

                if
                (
                    db().time().foundObject<uniformDimensionedVectorField>
                    (
                        cOfGdisplacement_ + "3Dof"
                    )
                )
                {
                    auto& rota =
                        db().time().lookupObjectRef<uniformDimensionedVectorField>
                        (
                            cOfGdisplacement_ + "3Dof"
                        );

                    rota.value()[0] = model_.cCofR(model_.bodyID(bodyIdCofG_))[0];
                    rota.value()[1] = model_.cCofR(model_.bodyID(bodyIdCofG_))[1];
                    rota.value()[2] += 0.5*t.deltaTValue()*(oldW + presentW);
                }
            }
            else
            {
                FatalErrorInFunction
                    << "CofGdisplacement is different to none." << endl
                    << "The model needs the entry body reference Id: bodyIdCofG."
                    << exit(FatalError);
            }
        }
    }

    forAll(forceFile_, bk)
    {
        /******************************************
        For two coordinates, the forces F and moemnts M (= r ^F) were expressed in one of them, you want to covert to another one, i.e., F' and M'
        you need some transformation information between two coordinates, such as rotation tensor A and translation vector P.
            For forces, F' = A & F,
            For moments, M' = (A & (r + P)) ^ (A & F) = A & (r ^ F + P ^ F) = A & (M + P ^ F)
        where rotation matrix and translation vector should be caculated according to the expression of the position of the new one in primitive one.
        but for fulfill the define spatialTransform::dual & spatialVector, the P was set to its opposite vector
        
        ******************************************/
        forceFile_[bk].calcForcesMoments();

        // for appendeges, express its forces and moments in its parent rigid body reference frame, if there is only one parent rigid
        //- for multi parent rigid body, it needs to be modified later.
        vector translation = model_.cCofR(model_.lambda()[bodyMeshes_[bk].bodyID_]); 
        // vector translation = model_.cCofR(bodyMeshes_[bk].bodyID_); 
        // Info<<nl<<"translation = "<<translation<<endl;

        scalar rotationA = 0.0;
        if
        (
            db().time().foundObject<uniformDimensionedVectorField>
            (
                cOfGdisplacement_ + "3Dof"
            )
        )
        {
            auto& rotationRef = db().time().lookupObject<uniformDimensionedVectorField>(cOfGdisplacement_ + "3Dof");
            rotationA = -1.0*rotationRef.value()[2]; // 
        }
        tensor rotationMatrix(cos(rotationA), -1*sin(rotationA), 0, sin(rotationA), cos(rotationA), 0, 0, 0, 1);
        spatialTransform transformDual(rotationMatrix, translation);
        
        spatialVector totalFM = spatialTransform::dual(transformDual) & spatialVector(forceFile_[bk].momentEff(), forceFile_[bk].forceEff());
        spatialVector pressFM = spatialTransform::dual(transformDual) & spatialVector(forceFile_[bk].pressureMomentEff(), forceFile_[bk].pressureForceEff());
        spatialVector viscoFM = spatialTransform::dual(transformDual) & spatialVector(forceFile_[bk].viscousMomentEff(), forceFile_[bk].viscousForceEff());

        scalar patchArea =  forceFile_[bk].patchArea();
        
        // spatialVector totalFMM = spatialTransform::dual(model_.X0(bodyMeshes_[bk].bodyID_)) & spatialVector(forceFile_[bk].momentEff(), forceFile_[bk].forceEff());
        if (Pstream::master() && model_.report())
        {
            Ostream& os = bodyDynamicsOs_[bk].file();
            bodyDynamicsOs_[bk].writeCurrentTime(os);
            os<<tab<<patchArea<<tab;
            for(label i =0; i< 6; i++)
                os<<totalFM[i]<<tab;
            for(label i =0; i< 6; i++)
                os<<pressFM[i]<<tab;
            for(label i =0; i< 6; i++)
                os<<viscoFM[i]<<tab;
            os<<endl;
        }
    }

    forAll(oumFile_, oi)
    {
        vector translation = model_.cCofR(model_.bodyID(bodyIdCofG_));
        scalar rotationA = 0.0;
        if
        (
            db().time().foundObject<uniformDimensionedVectorField>
            (
                cOfGdisplacement_ + "3Dof"
            )
        )
        {
            auto& rotationRef = db().time().lookupObject<uniformDimensionedVectorField>(cOfGdisplacement_ + "3Dof");
            rotationA = -1.0*rotationRef.value()[2]; // 
        }
        tensor rotationMatrix(cos(rotationA), -1*sin(rotationA), 0, sin(rotationA), cos(rotationA), 0, 0, 0, 1);
        spatialTransform transformDual(rotationMatrix, translation);
        spatialVector oumForce(Zero, Zero);
        if
        (
            mesh().time().foundObject<uniformDimensionedSymmTensorField>
            (
                oumName_[oi] + "sourceForce"
            )
        )
        
        {
            auto& fxx =
                mesh().time().lookupObjectRef<uniformDimensionedSymmTensorField>
                (
                    oumName_[oi] + "sourceForce"
                );
            for(label i =0; i < 6; i++)    
            oumForce[i] = fxx.value()[i];
        }
        spatialVector totalFM = spatialTransform::dual(transformDual) & oumForce;
        if (Pstream::master() && model_.report())
        {
            Ostream& os = oumFile_[oi].file();
            oumFile_[oi].writeCurrentTime(os);
            os<<tab;
            for(label i =0; i< 6; i++)
                os<<totalFM[i]<<tab;
            os<<endl;
        }
    }
   

    if (Pstream::master() && model_.report())
    {
    
        forAll(bodyMeshes_, bi)
        {
            model_.status(bodyMeshes_[bi].bodyID_);
        }

        scalarField qNew(bodyMeshes_.size()*6, Zero);
        scalarField qDotNew(bodyMeshes_.size()*6, Zero);
        scalarField qDdotNew(bodyMeshes_.size()*6, Zero);
        model_.writingState(model_.state().q(), model_.state().qDot(), model_.state().qDdot(), qNew, qDotNew, qDdotNew);

        forAll(bodyMotionsOs_, bj)
        {
            Ostream& os = bodyMotionsOs_[bj].file();
            bodyMotionsOs_[bj].writeCurrentTime(os);
            os<<tab;
            for(label i = bj*6; i< 6 + bj*6; i++)
                os<<qNew[i]<<tab;
            for(label i = bj*6; i< 6 + bj*6; i++)
                os<<qDotNew[i]<<tab;
            for(label i = bj*6; i< 6 + bj*6; i++)
                os<<qDdotNew[i]<<tab;
            os<<endl;
        }

        
    }
    if (t.value() >= rampTime_)
    {
        // Update the displacements
        if (bodyMeshes_.size() == 1)
        {
            pointDisplacement_.primitiveFieldRef() = model_.transformPoints
            (
                bodyMeshes_[0].bodyID_,
                bodyMeshes_[0].weight_,
                points0()
            ) - points0();
        }
        else
        {
            labelList bodyIDs(bodyMeshes_.size());
            List<const scalarField*> weights(bodyMeshes_.size());
            forAll(bodyIDs, bi)
            {
                bodyIDs[bi] = bodyMeshes_[bi].bodyID_;
                weights[bi] = &bodyMeshes_[bi].weight_;
            }

            pointDisplacement_.primitiveFieldRef() =
                model_.transformPoints(bodyIDs, weights, points0()) - points0();

        }
    }
    // Displacement has changed. Update boundary conditions
    pointConstraints::New
    (
        pointDisplacement_.mesh()
    ).constrainDisplacement(pointDisplacement_);
}


bool Foam::rigidBodyMeshMotion::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    // Force ASCII writing
    streamOpt.format(IOstream::ASCII);

    IOdictionary dict
    (
        IOobject
        (
            "rigidBodyMotionState",
            mesh().time().timeName(),
            "uniform",
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    IOdictionary maneuvers
    (
        IOobject
        (
            "maneuvers",
            mesh().time().timeName(),
            "uniform",
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    if(maneuveringOutput_.size())
    {
        forAll(maneuveringOutput_, bi)
        {
            maneuveringOutput_[bi].write(maneuvers);

        }
        maneuvers.regIOobject::writeObject(streamOpt, writeOnProc);
    }
    
    scalarField qNew(bodyMeshes_.size()*6, Zero);
    scalarField qDotNew(bodyMeshes_.size()*6, Zero);
    scalarField qDdotNew(bodyMeshes_.size()*6, Zero);
    model_.writingState(model_.state().q(), model_.state().qDot(), model_.state().qDdot(), qNew, qDotNew, qDdotNew);

    model_.state().write(dict, qNew, qDotNew, qDdotNew);
    return dict.regIOobject::writeObject(streamOpt, writeOnProc);
}


bool Foam::rigidBodyMeshMotion::read()
{
    if (displacementMotionSolver::read())
    {
        model_.read(coeffDict());

        return true;
    }

    return false;
}


// ************************************************************************* //
