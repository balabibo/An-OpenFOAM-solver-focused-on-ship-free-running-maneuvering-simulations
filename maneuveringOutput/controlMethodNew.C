//---------------------------------------------------------------------------*/
#include "controlMethod.H"
//---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * Selectors* * * * * * * * * * * * * * //


Foam::autoPtr<Foam::maneuvering::controlMethod> 
Foam::maneuvering::controlMethod::New
(
    const dictionary& dict
)
{
    const word motionType(dict.get<word>("type"));

    Info<< "Selecting maneuvering control type: "<< motionType <<"Control"<< endl;

    auto* ctorPtr = dictionaryConstructorTable(motionType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "controlMethod",
            motionType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<controlMethod>(ctorPtr(dict));
}  
