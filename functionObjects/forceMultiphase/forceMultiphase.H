/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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
    Foam::functionObjects::forceMultiphase

Group
    grpForcesFunctionObjects

Description
    Computes forces and moments over a given list of patches in a specific phase
    by integrating pressure and viscous forces and moments, and optionally resistance forces
    and moments from porous zones.

    Forces and moments are output in their total and constituent components:
    - total forces and moments
    - pressure contributions
    - viscous contributions
    - porous resistance contributions (optional)

    Forces and moments can be computed and output in:
    - the global Cartesian coordinate system (default)
    - a user-defined Cartesian coordinate system

    Operands:
    \table
      Operand       | Type           | Location
      input         | -              | -
      output file   | dat            | postProcessing/\<FO\>/\<time\>/\<file\>s
      output field  | volVectorField | \<time\>/\<outField\>s
    \endtable

    where \c \<file\>s:
    \verbatim
      force.dat        | Forces
      moment.dat       | Moments
    \endverbatim

    where \c \<outField\>s:
    \verbatim
      <namePrefix>:force   | Force field
      <namePrefix>:moment  | Moment field
    \endverbatim

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    <namePrefix>
    {
        // Mandatory entries
        type                forceMultiphase;
        libs                (forces);
        patches             (<wordRes>);

        // Optional entries
        directForceDensity  <bool>;
        porosity            <bool>;
        writeFields         <bool>;
        useNamePrefix       <bool>;

        // Conditional mandatory entries

            // if directForceDensity == true
            fD              <word>;


            // Cartesian coordinate system specification when
            // evaluating forces and moments, either of the below

            // Define the centre of rotation
            // with implicit directions e1=(1 0 0) and e3=(0 0 1)
            CofR                (0 0 0); // Centre of rotation

            // Define local coordinate system by origin + axes
            origin              (0 0 0);
            e1                  (1 0 0);
            e3                  (0 0 1); // (e1, e2) or (e2, e3) or (e3, e1)

            // General coordinate system specification (always cartesian)
            coordinateSystem
            {
                origin          (0 0 0);
                rotation
                {
                    type        axes;
                    e3          (0 0 1);
                    e1          (1 0 0); // (e1, e2) or (e2, e3) or (e3, e1)
                }
            }

        // Conditional optional entries

            // if directForceDensity == false
            p               <word>;
            U               <word>;
            phase           <word>;
            pRef            <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description               | Type | Reqd    | Deflt
      type       | Type name: forceMultiphase| word | yes     | -
      libs       | Library name: forces      | word | yes     | -
      patches    | Names of operand patches  | wordRes | yes  | -
      directForceDensity | Flag to directly supply force density <!--
                 -->                         | bool | no      | false
      porosity   | Flag to include porosity contributions | bool | no | false
      writeFields | Flag to write force and moment fields | bool | no | false
      useNamePrefix | Flag to include prefix for field names | bool | no | false
      coordinateSystem | Coordinate system specifier | dictionary | cndtnl | -
      CofR    | Centre of rotation          | vector | cndtnl   | -
      origin  | Origin of coordinate system | vector | cndtnl   | -
      e3      | e3 coordinate axis          | vector | cndtnl   | -
      e1      | e1 coordinate axis          | vector | cndtnl   | -
      fD      | Name of force density field | word   | cndtnl   | -
      p       | Name of pressure field      | word   | cndtnl   | p
      U       | Name of velocity field      | word   | cndtnl   | U
      phase   | Name of phase where forces are caculated | word | cndtnl | water
      pRef    | Value of reference pressure | scalar | cndtnl   | 0
    \endtable

    The inherited entries are elaborated in:
      - \link functionObject.H \endlink
      - \link writeFile.H \endlink
      - \link coordinateSystem.H \endlink

Note
  - For water-air multiphase flow case, if you want to caculate forces and moments
    only on patch immerged in water phase, you can set \c phase to \c water, or set
    \c phase to \c air to caculate the forces and moments in air phase only.
  - \c writeControl and \c writeInterval entries of function
    object do control when to output force and moment files and fields.
  - If a \c coordinateSystem entry exists, it is taken in favour of \c CofR.

SourceFiles
    forceMultiphase.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_forceMultiphase_H
#define Foam_functionObjects_forceMultiphase_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "coordinateSystem.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class forceMultiphase Declaration
\*---------------------------------------------------------------------------*/

class forceMultiphase
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected Data

        // Fields

            //- Sum of patch pressure forces
            vector sumPatchForcesP_;

            //- Sum of patch viscous forces
            vector sumPatchForcesV_;

            //- Sum of patch pressure moments
            vector sumPatchMomentsP_;

            //- Sum of patch viscous moments
            vector sumPatchMomentsV_;

            //- Sum of internal forces
            vector sumInternalForces_;

            //- Sum of internal moments
            vector sumInternalMoments_;


        // File streams

            //- File stream for forces
            autoPtr<OFstream> forceFilePtr_;

            //- File stream for moments
            autoPtr<OFstream> momentFilePtr_;


        // Read from dictionary

            //- Coordinate system used when evaluating forces and moments
            autoPtr<coordinateSystem> coordSysPtr_;

            //- Selected operand patches
            labelList patchIDs_;

            //- Reference density needed for incompressible calculations
            scalar rhoRef_;

            //- Reference pressure
            scalar pRef_;

            //- Name of pressure field
            word pName_;

            //- Name of velocity field
            word UName_;

            //- Name of density field
            word rhoName_;

            //- Name of force density field
            word fDName_;

            //- Flag to directly supply force density
            bool directForceDensity_;

            //- Flag to include porosity effects
            bool porosity_;

            //- Flag to write force and moment fields
            bool writeFields_;

            //- Flag of initialisation (internal)
            bool initialised_;

            //-Name of phase where forces and moments would be caculated
            word phaseName_;

            //-Area of patch where forces and moments would be caculated
            scalar patchArea_;


    // Protected Member Functions

        //- Set the co-ordinate system from dictionary and axes names
        void setCoordinateSystem
        (
            const dictionary& dict,
            const word& e3Name = word::null,
            const word& e1Name = word::null
        );

        //- Return access to the force field
        volVectorField& force();

        //- Return access to the moment field
        volVectorField& moment();

        //- Initialise containers and fields
        void initialise();

        //- Reset containers and fields
        void reset();


    // Evaluation

        //- Return the effective stress (viscous + turbulent) for patch
        tmp<symmTensorField> devRhoReff
        (
            const tensorField& gradUp,
            const label patchi
        ) const;

        //- Return dynamic viscosity field
        tmp<volScalarField> mu() const;

        //- Return rho if specified otherwise rhoRef
        tmp<volScalarField> rho() const;

        //- Return rho if specified otherwise rhoRef for patch
        tmp<scalarField> rho(const label patchi) const;

        //- Return rhoRef if the pressure field is
        //- dynamic (i.e. p/rho), otherwise return 1
        scalar rho(const volScalarField& p) const;

        //- Add patch contributions to force and moment fields
        void addToPatchFields
        (
            const label patchi,
            const vectorField& Md,
            const vectorField& fP,
            const vectorField& fV
        );

        //- Add cell contributions to force and
        //- moment fields, and include porosity effects
        void addToInternalField
        (
            const labelList& cellIDs,
            const vectorField& Md,
            const vectorField& f
        );


    // I-O

        //- Create the integrated-data files
        void createIntegratedDataFiles();

        //- Write header for an integrated-data file
        void writeIntegratedDataFileHeader
        (
            const word& header,
            OFstream& os
        ) const;

        //- Write integrated data to files
        void writeIntegratedDataFiles();

        //- Write integrated data to a file
        void writeIntegratedDataFile
        (
            const vector& pres,
            const vector& vis,
            const vector& internal,
            OFstream& os
        ) const;

        //- Write integrated data to stream
        void logIntegratedData
        (
            const string& descriptor,
            const vector& pres,
            const vector& vis,
            const vector& internal
        ) const;


public:

    //- Runtime type information
    TypeName("forceMultiphase");


    // Constructors

        //- Construct from Time and dictionary
        forceMultiphase
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool readFields = true
        );

        //- Construct from objectRegistry and dictionary
        forceMultiphase
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool readFields = true
        );

        //- No copy construct
        forceMultiphase(const forceMultiphase&) = delete;

        //- No copy assignment
        void operator=(const forceMultiphase&) = delete;


    //- Destructor
    virtual ~forceMultiphase() = default;


    // Member Functions

        //-Return the area of patch on specific phase
        scalar patchArea() const;
        
        //- Read the dictionary
        virtual bool read(const dictionary& dict);

        //- Calculate forces and moments
        virtual void calcForcesMoments();

        //- Return the total force
        virtual vector forceEff() const;

        //-Return the pressure force
        virtual vector pressureForceEff() const;

        //-Return the voscous force
        virtual vector viscousForceEff() const;

        //- Return the total moment
        virtual vector momentEff() const;

        //-Return the pressure force
        virtual vector pressureMomentEff() const;

        //-Return the voscous force
        virtual vector viscousMomentEff() const;

        //- Execute the function object
        virtual bool execute();

        //- Write to data files/fields and to streams
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
