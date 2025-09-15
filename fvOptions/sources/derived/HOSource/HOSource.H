/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 ENERCON GmbH
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

Class
    Foam::fv::HOSource

Group
    grpFvOptionsSources

Description
    Reference to G.R. Hough, D.E. Ordway The generalized actuator disk Develop. Theor. Appl. Mech., 2 (1964), p. 42
    
    An example to use 
        disk1
        {
            type         HOSource;
            refBody    hull; //let the disk moves along with the rigid-body
            actBody    disk1; //present disk name
            rotationDir  right; // or left
            KTCoeff (0.4873 -0.4399 0.0079 -0.0945 0 0);  // polynomial coefficients for KT-J: left to right, constant to 5th order
            KQCoeff (0.0694 -0.0561 0.0065 -0.0189 0 0);   // polynomial coefficients for KQ-J: left to right, constant to 5th order
            diskRadius     0.04637; 
            diskHub        0.18;
            diskDir         (-1 0 0);
            diskOrigin      (0 0 0);
            diskThickness    0.0064;
            diskRPS           5.071;
        
        }

SourceFiles
    HOSource.C
    oumSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef oumSource_H
#define oumSource_H

#include "fvCFD.H"
#include "fvOptions.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "geometricOneField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
//#include "cellSetOption.H"
#include "writeFile.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class HOSource Declaration
\*---------------------------------------------------------------------------*/

class HOSource
:
    public fv::option,
    public functionObjects::writeFile
{
protected:

     // specified in the dictionary

        //- Rotation direction, where -1 means left, 1 means right 
        label rotaDir_;

        //- Polynomial coefficients of the KT curve, 5 order max 
        scalarList KTCoef_;

        //- Polynomial coefficients of the KQ curve, 5 order max 
        scalarList KQCoef_;        
              
        //- The propeller radius
        scalar diskR_;

        //- propeller hub radius ratio realtive to propeller radius
        scalar diskH_;

        //- The vector of propeller axis
        vector diskDir_;

        //- The initial origin of propeller disk, without changes
        const vector oriDiskOri_;

        //- The origin of the propeller disk
        vector diskOri_;

        //- The abstract point P, which is equal to oriDiskOri_ + diskDir_
        //- It can be used to update the disk position 
        const vector p_;

        //- The thickness of the propeller disk
        scalar diskThick_;
               
        //- The initial rotation speed (rps)
        scalar diskRPS_;

        //- Name of reference rigid for dynamic mesh
        const word refBody_;

        //- Name of act momentum source for dynamic mesh
        const word actBody_;

        //- velocity of the rigid body which the source acts on 
        //-to acquire advance velocity for the source
        vector vcofR_;


    //- calculated
        
        // Set of cells to apply the source to
        DynamicList<label> cells_;

        // Set of cells to calculate advance coefficient J, considering wake
        DynamicList<label> wakeCells_;
        
        // Sum of cell volumes
        scalar diskV_;

        // Sum of wake cell volumes
        scalar wakeCellV_;

        //- external forces and moments due to current momentum source
        spatialVector fx_;
        
    // Protected Member Functions

        //- Output file header information
        //virtual void writeFileHeader(Ostream& os);


private:

    // Private Member Functions

         
         
        //- Compute the selected force computation method type
        template<class AlphaFieldType, class RhoFieldType>
        void calc
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<vector>& eqn,
            const labelList& cells
        );



public:

    //- Runtime type information
    TypeName("HOSource");


    // Constructors

        //- No default construct
        HOSource() = delete;

        //- Construct from components
        HOSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        HOSource(const HOSource&) = delete;

        //- No copy assignment
        void operator=(const HOSource&) = delete;


    //- Destructor
    virtual ~HOSource() = default;

    //Member function
    scalar KT(const scalar& J) const
    {
        return KTCoef_[0] + KTCoef_[1]*J + KTCoef_[2]*sqr(J) + KTCoef_[3]*pow3(J) + KTCoef_[4]*pow4(J) + KTCoef_[5]*pow5(J);
    }

    scalar KQ(const scalar& J) const
    {
        return KQCoef_[0] + KQCoef_[1]*J + KQCoef_[2]*sqr(J) + KQCoef_[3]*pow3(J) + KQCoef_[4]*pow4(J) + KQCoef_[5]*pow5(J);
    }    
     
    // Access

    //- const Normal disk direction
    const vector& diskDir() const
    {
        return diskDir_;
    }

    //- Normal disk direction
    vector& diskDir()
    {
        return diskDir_;
    }

    //- const Disk origin
    const vector& diskOri() const
    {
        return diskOri_;
    }

    //- Disk origin
    vector& diskOri()
    {
        return diskOri_;
    }

    //- const disk revolution per second
    const scalar& diskRPS() const
    {
        return diskRPS_;
    }

    //- disk revolution per second
    scalar& diskRPS()
    {
        return diskRPS_;
    }

    //- Disk volume
    const scalar& diskV() const
    {
        return diskV_;
    }

    //- reference rigid body
    const word refBody() const
    {
        return refBody_;
    }

    //- Disk volume
    const word actBody() const
    {
        return actBody_;
    }

    // Add explicit and implicit contributions

    //- Source term to momentum equation
    virtual void addSup
    (
        fvMatrix<vector>& eqn,
        const label fieldi
    );

    //- Source term to compressible momentum equation
    virtual void addSup
    (
        const volScalarField& rho,
        fvMatrix<vector>& eqn,
        const label fieldi
    );

    //- Source term to phase-compressible momentum equation
    virtual void addSup
    (
        const volScalarField& alpha,
        const volScalarField& rho,
        fvMatrix<vector>& eqn,
        const label fieldi
    );


    // IO

    //- Read dictionary
    virtual bool read(const dictionary& dict);


    // Update

    //- Update disk parameters
    void update();

    //- update source cells list
    void updateCells();

    //- update wake cells list
    // DynamicList<label> updateWakeCells();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HOSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
