/*---------------------------------------------------------------------------*\
Class
    Foam::maneuveringInput

Description
    Part of the maneuveringOutput.C， which is used to read sailing velocity or 
yaw angle of prescribed rigid body.

    using RTS function

\*---------------------------------------------------------------------------*/

#ifndef Foam_maneuveringInput_H
#define Foam_maneuveringInput_H

#include "fvCFD.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"
//---------------------------------------------------------------------------*/

namespace Foam
{
namespace maneuvering
{

class maneuveringInput
{

protected:
    
    const word inputName_;
    const word refBody_;
    const word actBody_;   

    // No copy construct and copy assignment
    maneuveringInput(const maneuveringInput&) = delete;
    void operator=(const maneuveringInput&) = delete;


public:

    TypeName("maneuveringInput");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            maneuveringInput,
            dictionary,
            (const dictionary& dict),
            (dict)
        ); 
        
    // Constructors

        //- Construct from dictionary
        maneuveringInput() = delete;
        maneuveringInput(const dictionary &dict);
        // virtual autoPtr<maneuveringInput> clone() const = 0;

    // Selectors

        //- Select constructed from the  ictionary 
        static autoPtr<maneuveringInput> New
        (
            const dictionary &dict
        );  

        virtual ~maneuveringInput() = default;

    // Member Functions


        const word inputName() const;
        const word refBody() const;
        const word actBody() const;

        virtual label inputTypeValue() const = 0;
        virtual label controlType() const = 0;

        // write 
        virtual void write(dictionary&) const;

};


class sailingInput 
: 
    public maneuveringInput
{
    public:
        
        TypeName ("sailingInput");

        sailingInput() = delete;
        sailingInput(const dictionary& dict);

        label inputTypeValue() const override;      
        label controlType() const override;

        void write(dictionary &) const override;
        
};


class yawInput 
: 
    public maneuveringInput
{
    public:
        TypeName ("yawInput");

        yawInput() = delete;
        yawInput(const dictionary &dict);

        label inputTypeValue() const override;
        label controlType() const override;

        void write(dictionary &) const override;

};

}
}
#endif
