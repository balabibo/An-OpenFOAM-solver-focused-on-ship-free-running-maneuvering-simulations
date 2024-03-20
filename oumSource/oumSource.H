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
    Applies sources on velocity (i.e. \c U) within a specified region
    to enable actuator disk models for aero/hydro thrust loading of
    horizontal axis turbines on surrounding flow field in terms of
    energy conversion processes.

    Available options for force computations:
    \verbatim
      Froude          | Froude's one-dimensional ideal actuator disk method
      variableScaling | Variable-scaling actuator disk method
    \endverbatim

    The expressions for \c Froude method:
    \f[
        T = 2 \rho A | \vec{U_o} \cdot \vec{n} |^2 a (1-a)
    \f]
    with
    \f[
        a = 1 - \frac{C_p}{C_T}
    \f]
    where
    \vartable
      T         | Thrust magnitude
      \rho      | Monitored incoming fluid density
      A         | Actuator disk planar surface area
      \vec{U_o} | Incoming velocity spatial-averaged on monitored region
      \vec{n}   | Surface-normal vector of the actuator disk pointing downstream
      a         | Axial induction factor
      C_p       | Power coefficient
      C_T       | Thrust coefficient
    \endvartable

    The expressions for \c variableScaling method:
    \f[
        T = 0.5 \rho A | \vec{U_m} \cdot \vec{n} |^2 C_T^*
    \f]
    with
    \f[
        C_T^* = C_T \left( \frac{ |\vec{U_{ref}}| }{ |\vec{U_m}| } \right)^2
    \f]
    where
    \vartable
      T         | Thrust magnitude
      \rho      | Monitored incoming fluid density
      A         | Actuator disk planar surface area
      \vec{U_m} | Incoming velocity spatial-averaged on actuator disk
      \vec{n}   | Surface-normal vector of the actuator disk pointing downstream
      \vec{U_{ref}} | Incoming velocity spatial-averaged on monitored region
      C_T       | Thrust coefficient
      C_T^*     | Calibrated thrust coefficient
    \endvartable

    Reference
    \verbatim
        Froude's one-dimensional ideal actuator disk (tags:F,BJSB):
            Froude, R. E. (1889).
            On the part played in propulsion by differences of fluid pressure.
            Trans. Inst. Naval Architects, 30, 390.

            Burton, T., Jenkins, N., Sharpe, D., & Bossanyi, E. (2011).
            Wind energy handbook, 2nd edition.
            Chichester, the United Kingdom. John Wiley & Sons.
            Print ISBN:978-0-470-69975-1
            DOI:10.1002/9781119992714

        Variable-scaling method (tags:LSRMKT,LS):
            van der Laan, M. P., Sørensen, N. N., Réthoré, P. E.,
            Mann, J., Kelly, M. C., & Troldborg, N. (2015).
            The k‐ε‐fP model applied to double wind turbine
            wakes using different actuator disk force methods.
            Wind Energy, 18(12), 2223-2240.
            DOI:10.1002/we.1816

            van der Laan, M. P., Sørensen, N. N., Réthoré, P. E., Mann, J.,
            Kelly, M. C., Troldborg, N., Hansen, K. S., & Murcia, J. P. (2015).
            The k‐ε‐fP model applied to wind farms.
            Wind Energy, 18(12), 2065-2084.
            DOI:10.1002/we.1804
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
    oumSource1
    {
        // Mandatory entries (unmodifiable)
        type            oumSource;

        // Mandatory entries (runtime modifiable)
        diskArea        40.0;
        diskDir         (1 0 0);
        Cp              <Function1>;
        Ct              <Function1>;

        // Conditional optional entries (unmodifiable)
        monitorMethod   <points|cellSet>;
        monitorCoeffs
        {
            // Option-1
            points
            (
                (p1x p1y p1z)
                (p2x p2y p2z)
                ...
            );

            // Option-2
            cellSet     <monitorCellSet>;
        }

        // Optional entries (unmodifiable)
        variant         <forceMethod>;

        // Optional entries (runtime modifiable)
        sink            true;
        writeFileStart  0;
        writeFileEnd    100;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Reqd     | Dflt
      type       | Type name: oumSource      | word | yes      | -
      diskArea   | Actuator disk planar surface area   | scalar | yes    | -
      diskDir    | Surface-normal vector of the actuator disk <!--
               --> pointing downstream                 | vector    | yes | -
      Cp         | Power coefficient                   | Function1 | yes | -
      Ct         | Thrust coefficient                  | Function1 | yes | -
      monitorMethod | Type of incoming velocity monitoring method   <!--
               -->    - see below                      | word | no  | points
      variant    | Type of the force computation method - see below <!--
               -->                                     | word | no  | Froude
      sink       | Flag for body forces to act as a source (false)  <!--
               -->  or a sink (true)                   | bool   | no  | true
      writeFileStart | Start time for file output      | scalar | no  | 0
      writeFileEnd   | End time for file output        | scalar | no  | VGREAT
    \endtable

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink
      - \link writeFile.H \endlink
      - \link Function1.H \endlink

    Options for the \c monitorMethod entry:
    \verbatim
      points  | Monitor incoming velocity field at a given set of points
      cellSet | Monitor incoming velocity field at a given cellSet
    \endverbatim

    Options for the \c variant entry:
    \verbatim
      Froude          | Froude's one-dimensional ideal actuator disk method
      variableScaling | Variable-scaling actuator disk method
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
//#include "cellSetOption.H"
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

        //- Center of rotation of the rigid body which the source acts on
        vector cofR_;


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
