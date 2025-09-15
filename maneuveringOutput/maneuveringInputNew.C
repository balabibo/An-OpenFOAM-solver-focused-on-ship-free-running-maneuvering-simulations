/*---------------------------------------------------------------------------*\
Class
    Foam::maneuveringInput

Description
    Part of the maneuveringOutput.C， which is used to read sailing velocity or 
yaw angle of prescribed rigid body.

using RTS function

\*---------------------------------------------------------------------------*/

#include "maneuveringInput.H"


// * * * * * * * * * * * * * * * * Selectors* * * * * * * * * * * * * * //


Foam::autoPtr<Foam::maneuvering::maneuveringInput> 
Foam::maneuvering::maneuveringInput::New
(
    const dictionary &dict
)
{
    const word inputType(dict.get<word>("type"));
    word motionType;
    
    if(inputType == "sailing")
    {
        motionType = "sailingInput";
    }
    
    else
    {
        motionType = "yawInput";
    }

    Info<< "Selecting maneuvering input type: "<< motionType << endl;

    auto* ctorPtr = dictionaryConstructorTable(motionType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "maneuveringInput",
            motionType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<maneuveringInput>(ctorPtr(dict));
}  
