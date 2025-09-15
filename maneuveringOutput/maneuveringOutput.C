#include "maneuveringOutput.H"


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::maneuveringOutput::maneuveringOutput(const polyMesh &mesh, const dictionary &dict)
: 
    mesh_(mesh),
    maneuveringInput_(Foam::maneuvering::maneuveringInput::New(dict)),
    controlMethod_(Foam::maneuvering::controlMethod::New(dict)),
    timeIndex_(mesh.time().timeIndex())
{}


// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * *//

scalar Foam::maneuveringOutput::output(const vector2D& input)
{
    // Get time data
    const scalar deltaT = mesh_.time().deltaTValue();
    const scalar t = mesh_.time().timeOutputValue();
    if(t <  controlMethod_.get()->cStartTime() || t > controlMethod_.get()->cEndTime())
    {
       return 0.0; 
    }

    // Update the old-time quantities
    if (timeIndex_ != mesh_.time().timeIndex())
    {
        timeIndex_ = mesh_.time().timeIndex();
    }
    //- acquire input value
    const scalar inputValue = input[maneuveringInput_.get()->inputTypeValue()];
        
    const scalar outputSignal = controlMethod_.get()->calculate(inputValue, deltaT);
    

    //Info << "maneuveringOutput: targetValue = " << targetValue_->value(t) << endl;
    Info << "maneuveringOutput: input value = " << inputValue << endl;
    //Info << "maneuveringOutput: error = " << targetValue_->value(t) - maneuveringInputValue  << endl;
    Info << "maneuveringOutput: outputSignal = " << outputSignal << endl;

    return outputSignal;
}

void Foam::maneuveringOutput::write(Ostream& os) const
{


}

void Foam::maneuveringOutput::write(dictionary& dict) const
{

    dictionary& enDict = dict.subDictOrAdd(controlMethod_.get()->controlName());
    maneuveringInput_.get()->write(enDict);
    controlMethod_.get()->write(enDict);

}


const Foam::maneuvering::maneuveringInput& maneuveringOutput::mInput() const
{
    return maneuveringInput_();
}







