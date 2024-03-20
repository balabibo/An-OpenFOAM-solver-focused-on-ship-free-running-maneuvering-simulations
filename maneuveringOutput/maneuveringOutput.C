#include "maneuveringOutput.H"
#include "IOstreams.H"


// * * * * * * * * * * * * Base maneuveringOutput  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

maneuveringOutput::maneuveringOutput(const polyMesh &mesh, const dictionary &dict)
: 
    mesh_(mesh),
    maneuveringInput_(maneuveringInput::create(dict)),
    controlMethod_(controlMethod::create(dict)),
    timeIndex_(mesh.time().timeIndex())
{}


// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * *//

scalar maneuveringOutput::output(const vector2D& input)
{
    // Get time data
    const scalar deltaT = mesh_.time().deltaTValue();
    const scalar t = mesh_.time().timeOutputValue();
    if(t <  controlMethod_->cStartTime() || t > controlMethod_->cEndTime())
    {
       return 0.0; 
    }

    // Update the old-time quantities
    if (timeIndex_ != mesh_.time().timeIndex())
    {
        timeIndex_ = mesh_.time().timeIndex();
    }
    //- acquire input value
    const scalar inputValue = input[maneuveringInput_->inputTypeValue()];
        
    const scalar outputSignal = controlMethod_->calculate(inputValue, deltaT);

    //Info << "maneuveringOutput: targetValue = " << targetValue_->value(t) << endl;
    Info << "maneuveringOutput: input value = " << inputValue << endl;
    //Info << "maneuveringOutput: error = " << targetValue_->value(t) - maneuveringInputValue  << endl;
    Info << "maneuveringOutput: outputSignal = " << outputSignal << endl;

    return outputSignal;
}

void maneuveringOutput::write(Ostream& os) const
{


}

void maneuveringOutput::write(dictionary& dict) const
{

    dictionary& enDict =dict.subDictOrAdd(controlMethod_->controlName());
    maneuveringInput_->write(enDict);
    controlMethod_->write(enDict);

}


const std::shared_ptr<maneuveringInput> maneuveringOutput::mInput()
{
    return maneuveringInput_;
}







