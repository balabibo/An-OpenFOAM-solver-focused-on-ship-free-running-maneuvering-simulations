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
    Foam::fv::oumSource

Group
    grpFvOptionsSources

Description
    A momentum source based on blade element momentum theory, you can find
    more details on the reference.
   
    Reference
    \verbatim
        
            Tokgoz Emel. (2015).
            A CFD study on the propeller-hull interaction flow in waves 
            using body-force distribution model
            DOI:10.18910/53979

    \endverbatim

    Sources applied to:
    \verbatim
      U         | Velocity    [m/s]
    \endverbatim

    Required fields:
    \verbatim
      U         | Velocity    [m/s]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    disk1
    {

      type            oumSource;
      refBody         hull;  // reference rigid body when you use maneuvering controller
      actBody         disk1; // act momentum souce when you use maneuvering controller
      chordCoeff (6.2575 -36.837 125.59 -226.05 201.38 -70.321);  // coefficients of a chord fitting function, form the constant term to the fifth order term.
      pitchCoeff (21.437 -108.69 295.87 -430.02 312.63 -89.56);   // coefficients of a pitch fitting function, form the constant term to the fifth order term.
      diskRadius    0.0465; // radius of propeller
      diskHub        0.186; // relative hub, which is equal to hub(unit: m)/diskRadius(unit: m)
      diskDir         (-1 0 0); // the normal direction of the propeller, and the revolution direction is decided according to right-hand law.
      diskOrigin     (2.65275 0 -0.078651); // the center position of the propeller
      diskThickness    0.01897182; // the thickness of the propeller
      bladeNumber             5; // number of blades
      diskRPS            19.31970308; // initial revolution per second
   
    }
    \endverbatim
    
SourceFiles
    oumSource.C
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
#include "writeFile.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class oumSource Declaration
\*---------------------------------------------------------------------------*/

class oumSource
:
    public fv::option,
    public functionObjects::writeFile
{
protected:

     // specified in the dictionary

        //- Rotation direction, where -1 means left, 1 means right 
        label rotaDir_;
        
        //- Polynomial coefficients of the chord length fit curve, 5 order max 
        scalarList chordCoef_;

        //- Polynomial coefficients of the pitch fit curve, 5 order max 
        scalarList pitchCoef_;
              
        //- The propeller radius
        scalar diskR_;

        //- propeller hub radius realtive to propeller radius
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
        
        //- The number of blades
        label bladeN_;
        
        //- The initial rotation speed (rps)
        scalar diskRPS_;

        //- Name of reference rigid for dynamic mesh
        const word refBody_;

        //- Name of act momentum source for dynamic mesh
        const word actBody_;

    //- calculated
        
        // Set of cells to apply the source to
        DynamicList<label> cells_;
        
        // Sum of cell volumes
        scalar diskV_;

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
    TypeName("oumSource");


    // Constructors

        //- No default construct
        oumSource() = delete;

        //- Construct from components
        oumSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        oumSource(const oumSource&) = delete;

        //- No copy assignment
        void operator=(const oumSource&) = delete;


    //- Destructor
    virtual ~oumSource() = default;


    // Member Functions
    
     template <typename T> T chordLength(T r)
    {
      return chordCoef_[5]*pow(r,5)+chordCoef_[4]*pow(r,4)+chordCoef_[3]*pow(r,3)+chordCoef_[2]*pow(r,2)+chordCoef_[1]*r+chordCoef_[0];
    }
    
    template <typename T> T pitchValue(T r)
    {
      return pitchCoef_[5]*pow(r,5)+pitchCoef_[4]*pow(r,4)+pitchCoef_[3]*pow(r,3)+pitchCoef_[2]*pow(r,2)+pitchCoef_[1]*r+pitchCoef_[0];
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

    //- update cells list
    DynamicList<label> updateCells();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "oumSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
