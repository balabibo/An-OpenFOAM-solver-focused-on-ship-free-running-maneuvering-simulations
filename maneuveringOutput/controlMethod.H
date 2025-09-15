/*---------------------------------------------------------------------------*\
Class
    Foam::controlMethod

Description
    controlMethod is a base class for implementing a control algotirhm
    used by maneuveringOutout.C.

    using RTS 

\*---------------------------------------------------------------------------*/

#ifndef controlMethod_H
#define controlMethod_H

#include "fvCFD.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"
//---------------------------------------------------------------------------*/

namespace Foam
{
namespace maneuvering
{

class controlMethod
{

protected:
    
    // Name of control Type
    word controlType_;

    // Start time of controller
    scalar cStartTime_;
    
    // End time of controller
    scalar cEndTime_;

    //Name of control module
    word controlName_;   

public:

    TypeName("controlMethod");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            controlMethod,
            dictionary,
            (const dictionary& dict),
            (dict)
        ); 
        
    // Constructors

        //- Construct from dictionary
        controlMethod(const dictionary &dict);
        controlMethod() = delete;

    // Selectors

        //- Select constructed from the  ictionary 
        static autoPtr<controlMethod> New
        (
            const dictionary &dict
        );  

        virtual ~controlMethod() = default;
    // Member Functions

        //- Calculate output signal from current process value and current time delta
        virtual scalar calculate(scalar, const scalar ) const = 0;

        //- Write to stream
        virtual void write(Ostream &) const;

        //- Write to dictionary
        virtual void write(dictionary &) const;    
        
        //- refer to start Time
        scalar cStartTime() const;
        
        //- refer to end Time
        scalar cEndTime() const;
        
        //- refer to name of control module
        const word controlName() const; 
};

class turningControl 
: 
    public controlMethod
{
    public:

        TypeName ("turning");

        turningControl() = delete;

        turningControl(const dictionary& dict);
        
        // currentYaw refers to current Yaw angle

        scalar calculate(scalar currentYaw, const scalar deltaT) const;

        void write(Ostream &) const override;

        void write(dictionary &) const override;

    private:

        // target yaw angle
        mutable scalar cTarget_;
        
        // maximum controller value 
        mutable scalar cMax_;
        
        // controller rate  
        mutable scalar cRate_;
        
        // Preserve old output signal for histeresis correction
        mutable scalar outputSignal_;

};

class zigzagControl 
: 
    public controlMethod
{
    public:

        TypeName ("zigzag");

        zigzagControl() = delete;

        zigzagControl(const dictionary &);
        
        // currentYaw refers to current yaw angle.
        scalar calculate(scalar currentYaw, scalar deltaT) const;
        
        void write(Ostream &) const override;

        void write(dictionary &) const override;

    private:
        
        // target yaw angle
        mutable scalar cTarget_;
        
        // maximum controller value 
        mutable scalar cMax_;
        
        // controller rate  
        mutable scalar cRate_;
        
        // output rudder angle
        mutable scalar outputSignal_;
        
        // restore the yaw angle of the previous time-step.
        mutable scalar oldYaw_; 
    
};

class sailingControl 
: 
    public controlMethod
{
    public:

        TypeName ("sailing");

        sailingControl() = delete;

        sailingControl(const dictionary &);
    
        // currentV refers to current sailing velocity
        scalar calculate(scalar currentV, scalar deltaT) const;
        
        void write(Ostream &) const override;

        void write(dictionary &) const override;


    private:
        mutable scalar P_;                     // Proportional control term
        mutable scalar I_;                     // Integral control term
        mutable scalar D_;                     // Differential control term
        mutable scalar cTarget_;                // target sailing velocity
        mutable scalar outputMax_;              // Max output signal (saturation)
        mutable scalar outputMin_;              // Min output signal
        mutable scalar errorMax_;               // Constains received absolute error value
        mutable scalar integralErrorMax_;       // Constains built-up integral error
        mutable scalar oldError_;               // Old error
        mutable scalar errorIntegral_;          // Error integral w.r.t. time
        mutable scalar outputSignal_;           // output RPS
};

class coursekeepingControl 
: 
    public controlMethod
{
    public:

        TypeName ("coursekeeping");

        coursekeepingControl() = delete;

        coursekeepingControl(const dictionary &);
    
        // currentV refers to current sailing velocity
        scalar calculate(scalar currentYaw, scalar deltaT) const;
        
        void write(Ostream &) const override;

        void write(dictionary &) const override;


    private:
        mutable scalar P_;                     // Proportional control term
        mutable scalar I_;                     // Integral control term
        mutable scalar D_;                     // Differential control term
        mutable scalar cTarget_;                // target heading angle
        mutable scalar maxRate_;                // maximum rudder rate
        mutable scalar maxAngle_;               // maximum rudder angle

        mutable scalar errorIntegral_;          // Error integral w.r.t. time
        mutable scalar oldError_;                 // old error
        mutable scalar oldRudderRate_;            // old rudder Rate
        mutable scalar oldRudderAngle_;           // old rudder angle
};

}
}
#endif
