#include "maneuveringInput.H"

const Foam::Enum
<
maneuveringInput::inputType
>
maneuveringInput::inputTypeNames
({
        {inputType::sailing, "sailing"},
        {inputType::turning, "turning"},
        {inputType::zigzag, "zigzag"},
        {inputType::coursekeeping, "coursekeeping"},
});

std::shared_ptr<maneuveringInput>
maneuveringInput::create(const dictionary &dict)
{
    const inputType type = inputTypeNames.get(dict.get<word>("type"));

    switch (type)
    {
    case sailing:
        return std::make_shared<sailingInput>(dict);
    case turning:
        return std::make_shared<yawInput>(dict);
    case zigzag:
        return std::make_shared<yawInput>(dict);     
    case coursekeeping:
        return std::make_shared<yawInput>(dict);          
    default:
        FatalIOErrorInFunction(dict)
            << "    Unknown control method " << type
            << exit(FatalIOError);
        return nullptr;
    }
}



// * * * * * * * * * * * * Base maneuveringInput  * * * * * * * * * * * * //

maneuveringInput::maneuveringInput(const dictionary &dict)
:
    inputName_(dict.get<word>("type")),
    refBody_(dict.get<word>("refBody")),
    actBody_(dict.get<word>("actBody"))
{
    Info<<nl<<"reference rigid body："<<refBody_<<nl
        <<"acting rigid body："<<actBody_<<endl;
}



const word maneuveringInput::inputName() const
{
    return inputName_;
}

const word maneuveringInput::refBody() const
{
    return refBody_;
}

const word maneuveringInput::actBody() const
{
    return actBody_;   
}

void maneuveringInput::write(dictionary& dict) const
{
    dict.add("type", inputName_);
    dict.add("refBody", refBody_);
    dict.add("actBody", actBody_);
}

// * * * * * * * * * * * * sailingInput  * * * * * * * * * * * * //
sailingInput::sailingInput(const dictionary &dict)
:
    maneuveringInput(dict)
{
  
}

label sailingInput::inputTypeValue() const
{ 
   return 0;
}

label sailingInput::controlType() const
{
   return 0;
}

void sailingInput::write(dictionary& dict) const
{
    maneuveringInput::write(dict);
}

// * * * * * * * * * * * * yawInput  * * * * * * * * * * * * //
yawInput::yawInput(const dictionary &dict)
:
    maneuveringInput(dict)
{
    
}

label yawInput::inputTypeValue() const
{ 
   return 1;
}

label yawInput::controlType() const
{
   return 2;
}

void yawInput::write(dictionary& dict) const
{
    maneuveringInput::write(dict);
}


