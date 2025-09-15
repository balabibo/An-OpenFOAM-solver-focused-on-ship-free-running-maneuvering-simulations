/*---------------------------------------------------------------------------*\

\*---------------------------------------------------------------------------*/

#include "maneuveringInput.H"
//---------------------------------------------------------------------------*/


namespace Foam
{
namespace maneuvering
{
    defineTypeNameAndDebug(maneuveringInput, 0);
    defineRunTimeSelectionTable(maneuveringInput, dictionary);

    defineTypeNameAndDebug(sailingInput, 0);
    addToRunTimeSelectionTable
    (
        maneuveringInput,
        sailingInput,
        dictionary
    );

    defineTypeNameAndDebug(yawInput, 0);
    addToRunTimeSelectionTable
    (
        maneuveringInput,
        yawInput,
        dictionary
    );

}
}

// * * * * * * * * * * * * * * * * base * * * * * * * * * * * * * * //

Foam::maneuvering::maneuveringInput::maneuveringInput(const dictionary &dict)
:
    inputName_(dict.get<word>("type")),
    refBody_(dict.get<word>("refBody")),
    actBody_(dict.get<word>("actBody"))
{
    Info<<nl<<"reference rigid body: "<<refBody_<<nl
        <<"acting rigid body: "<<actBody_<<endl;
}

const word Foam::maneuvering::maneuveringInput::inputName() const
{
    return inputName_;
}

const word Foam::maneuvering::maneuveringInput::refBody() const
{
    return refBody_;
}

const word Foam::maneuvering::maneuveringInput::actBody() const
{
    return actBody_;   
}

void Foam::maneuvering::maneuveringInput::write(dictionary& dict) const
{
    dict.add("type", inputName_);
    dict.add("refBody", refBody_);
    dict.add("actBody", actBody_);
}

// * * * * * * * * * * * * sailingInput  * * * * * * * * * * * * //
Foam::maneuvering::sailingInput::sailingInput(const dictionary &dict)
:
    maneuveringInput(dict)
{
  
}

label Foam::maneuvering::sailingInput::inputTypeValue() const
{ 
   return 0;
}

label Foam::maneuvering::sailingInput::controlType() const
{
   return 0;
}

void Foam::maneuvering::sailingInput::write(dictionary& dict) const
{
    maneuveringInput::write(dict);
}

// * * * * * * * * * * * * yawInput  * * * * * * * * * * * * //
Foam::maneuvering::yawInput::yawInput(const dictionary &dict)
:
    maneuveringInput(dict)
{
    
}

label Foam::maneuvering::yawInput::inputTypeValue() const
{ 
   return 1;
}

label Foam::maneuvering::yawInput::controlType() const
{
   return 2;
}

void Foam::maneuvering::yawInput::write(dictionary& dict) const
{
    maneuveringInput::write(dict);
}


