#include "controlMethod.H"

// * * * * * * * * * * * * Utility function  * * * * * * * * * * * * //



// * * * * * * * * * * * * Factory  * * * * * * * * * * * * //
const Foam::Enum
<
controlMethod::controlType
>
controlMethod::controlTypeNames
({
        {controlType::sailing, "sailing"},
        {controlType::turning, "turning"},
        {controlType::zigzag, "zigzag"},
        {controlType::coursekeeping, "coursekeeping"},
});

std::shared_ptr<controlMethod>
controlMethod::create(const dictionary &dict)
{
    const controlType type = controlTypeNames.get(dict.get<word>("type"));

    switch (type)
    {
    case sailing:
        return std::make_shared<sailingControl>(dict);
    case turning:
        return std::make_shared<turningControl>(dict);
    case zigzag:
        return std::make_shared<zigzagControl>(dict);
    case coursekeeping:
        return std::make_shared<coursekeepingControl>(dict);           
    default:
        FatalIOErrorInFunction(dict)
            << "    Unknown control method " << type
            << exit(FatalIOError);
        return nullptr;
    }
}

// * * * * * * * * * * * * Constructor  * * * * * * * * * * * * //
controlMethod::controlMethod(const dictionary &dict)
: 
  controlType_
  (
    controlTypeNames.get(dict.get<word>("type"))
  ),
  cStartTime_(dict.getOrDefault<scalar>("controllerStartTime", 0.)),
  cEndTime_(dict.getOrDefault<scalar>("controllerEndTime", 10000.)),
  controlName_(dict.dictName())
{}

void controlMethod::write(Ostream &os) const
{
    os.writeEntry("maneuveringMode", controlTypeNames.get(controlType_));
}

// * * * * * * * * * * * * member function  * * * * * * * * * * * * //
void controlMethod::write(dictionary& dict) const
{
    dict.add("controllerStartTime", cStartTime_);
    dict.add("controllerEndTime", cEndTime_);
}

scalar controlMethod::cStartTime() const
{
    return cStartTime_;
}

scalar controlMethod::cEndTime() const
{
    return cEndTime_;
}

enum controlMethod::controlType controlMethod::cType() const
{
    return controlType_;
}

const word controlMethod::controlName() const
{
    return controlName_;
}

// * * * * * * * * * * * * turning Control  * * * * * * * * * * * * //
turningControl::turningControl(const dictionary &dict)
: 
    controlMethod(dict),
    cTarget_(dict.getOrDefault<scalar>("controllerYawAngle", 270.)),
    cMax_(dict.getOrDefault<scalar>("controllerRudder", 35.)),
    cRate_(dict.getOrDefault<scalar>("controllerRate", 5.)),
    outputSignal_(dict.getOrDefault<scalar>("outputSignal", 0.))
{}

scalar turningControl::calculate(scalar currentYaw, const scalar deltaT)
{
    Info<<nl<<"turningControl is running!!"<<nl;
    if(abs(outputSignal_) >= abs(cMax_)*M_PI/180)
    {
       outputSignal_ = cMax_*M_PI/180;
       //return outputSignal_;
       return 0.0;
    }
    outputSignal_ += cRate_*deltaT*M_PI/180;
    //return outputSignal_;
    return cRate_*M_PI/180;
}

void turningControl::write(Ostream &os) const
{
    controlMethod::write(os);
    os.beginBlock("parameters");
    os.writeEntryIfDifferent("controllerTarget", 0.,  cTarget_);
    os.endBlock();
}

void turningControl::write(dictionary& dict) const
{
    controlMethod::write(dict);
    dict.add("controllerYawAngle", cTarget_);
    dict.add("controllerRudder", cMax_);
    dict.add("controllerRate", cRate_);
    dict.add("outputSignal", outputSignal_);
    
}

// * * * * * * * * * * * * zigzag Control  * * * * * * * * * * * * //
zigzagControl::zigzagControl(const dictionary &dict)
: 
    controlMethod(dict),
    cTarget_(dict.getOrDefault<scalar>("controllerYawAngle", 20.)),
    cMax_(dict.getOrDefault<scalar>("controllerRudder", 20.)),
    cRate_(dict.getOrDefault<scalar>("controllerRate", 5.)),
    outputSignal_(dict.getOrDefault<scalar>("outputSignal", 0.)),
    oldYaw_(dict.getOrDefault<scalar>("oldYaw", 0.))
{
 
}

scalar zigzagControl::calculate(scalar currentYaw, const scalar deltaT)
{
    Info<<nl<<"zigzagControl is running!!"<<nl;
    if(oldYaw_ < abs(cTarget_) && currentYaw >= abs(cTarget_))
    {
       cRate_ = -1*abs(cRate_);
    }
       
    else if(oldYaw_ > -1*abs(cTarget_) && currentYaw <= -1*abs(cTarget_))
    {
       cRate_ = abs(cRate_);      
    } 
       
    outputSignal_ += cRate_*deltaT;
    
    if(outputSignal_ >= abs(cMax_))
    {
       outputSignal_ = abs(cMax_);
       oldYaw_ = currentYaw;
       return 0.0;    
    }
    
    else if(outputSignal_ <= -1*abs(cMax_))
    {
       outputSignal_ = -1*abs(cMax_);
       oldYaw_ = currentYaw;
       return 0.0;        
    }
          
    oldYaw_ = currentYaw;
    
    //return outputSignal_;
    return cRate_;
}

void zigzagControl::write(Ostream &os) const
{
    controlMethod::write(os);
    os.beginBlock("parameters");
    os.writeEntryIfDifferent("controllerTarget", 0.,  cTarget_);
    os.endBlock();
}

void zigzagControl::write(dictionary& dict) const
{
    controlMethod::write(dict);
    dict.add("controllerYawAngle", cTarget_);
    dict.add("controllerRudder", cMax_);
    dict.add("controllerRate", cRate_);
    dict.add("outputSignal", outputSignal_);
    dict.add("oldYaw", oldYaw_);
    
}

// * * * * * * * * * * * * sailing Control  * * * * * * * * * * * * //
sailingControl::sailingControl(const dictionary &dict)
: 
    controlMethod(dict),
    P_(dict.getOrDefault<scalar>("controllerP", 1.)),
    I_(dict.getOrDefault<scalar>("controllerI", 1.)),
    D_(dict.getOrDefault<scalar>("controllerD", 0.)),
    cTarget_(dict.getOrDefault<scalar>("controllerTarget", 1.)),
    outputMax_(dict.getOrDefault<scalar>("controllerMax", 100.)),
    outputMin_(dict.getOrDefault<scalar>("controllerMin", 1.)),
    errorMax_(16.),
    integralErrorMax_(VGREAT),
    oldError_(dict.getOrDefault<scalar>("oldError", 0.)),
    errorIntegral_(dict.getOrDefault<scalar>("errorIntegral", 0.)),
    outputSignal_(dict.getOrDefault<scalar>("outputSignal", 0.))
{
    Info<<nl<<"********************************"<<nl
    <<"oldError: "<<oldError_<<nl
    <<"errorIntegral: "<<errorIntegral_<<nl
    <<"outputSignal: "<<outputSignal_<<endl;
}

scalar sailingControl::calculate(scalar currentV, scalar deltaT)
{
    Info<<nl<<"sailingControl is running!!"<<nl;
    Info<<nl<<"errorIntegral: "<<errorIntegral_<<endl;
    scalar error = fabs(cTarget_) - fabs(currentV);
    error = max(min(error, errorMax_), -errorMax_);  // Constain error according to specified errorMax
    errorIntegral_ += error * deltaT;
    errorIntegral_ = max(min(errorIntegral_, integralErrorMax_), -integralErrorMax_);
    const scalar errorDifferential = (error - oldError_);
    oldError_ = error;

    // Calculate increased output RPS value

    const scalar increasedOutputSignal = P_*error + I_*errorIntegral_ + D_*errorDifferential;
    outputSignal_ = increasedOutputSignal;
    
    // Return result within defined regulator saturation: outputMax_ and outputMin_
    return max(min(outputSignal_, outputMax_), outputMin_);
}

void sailingControl::write(Ostream &os) const
{
    controlMethod::write(os);
    os.beginBlock("parameters");
    os.writeEntry("Kp", P_);
    os.writeEntry("Ti", I_);
    os.writeEntry("Td", D_);
    os.writeEntryIfDifferent("outputMax", 1., outputMax_);
    os.writeEntryIfDifferent("outputMin", 0., outputMin_);
    os.writeEntryIfDifferent("errMax", VGREAT, errorMax_);
    os.writeEntryIfDifferent("errIntegMax", VGREAT, integralErrorMax_);
    os.endBlock();
}

void sailingControl::write(dictionary& dict) const
{
    controlMethod::write(dict);
    dict.add("controllerP", P_);
    dict.add("controllerI", I_);
    dict.add("controllerD", D_);
    dict.add("controllerTarget", cTarget_);
    dict.add("controllerMax", outputMax_);
    dict.add("controllerMin", outputMin_);
    dict.add("oldError", oldError_);
    dict.add("errorIntegral", errorIntegral_);
    dict.add("outputSignal", outputSignal_);
}

// * * * * * * * * * * * * coursekeeping Control  * * * * * * * * * * * * //
coursekeepingControl::coursekeepingControl(const dictionary &dict)
: 
    controlMethod(dict),
    P_(dict.getOrDefault<scalar>("controllerP", 1.)),
    I_(dict.getOrDefault<scalar>("controllerI", 1.)),
    D_(dict.getOrDefault<scalar>("controllerD", 0.)),
    cTarget_(dict.getOrDefault<scalar>("controllerTarget", 0.)),
    cRate_(dict.getOrDefault<scalar>("controllerRate", 5.)),
    outputMax_(dict.getOrDefault<scalar>("controllerMax", 35.)),// maximum rudder rate
    outputMin_(dict.getOrDefault<scalar>("controllerMin", -35.)),//minimum rudder rate
    errorMax_(16.),
    integralErrorMax_(VGREAT),
    oldError_(dict.getOrDefault<scalar>("oldError", 0.)),
    errorIntegral_(dict.getOrDefault<scalar>("errorIntegral", 0.)),
    outputSignal_(dict.getOrDefault<scalar>("outputSignal", 0.))
{}

scalar coursekeepingControl::calculate(scalar currentYaw, scalar deltaT)
{
    Info<<nl<<"coursekeepingControl is running!!"<<nl;
    scalar error = cTarget_ - currentYaw;
    error = max(min(error, errorMax_), -errorMax_);  // Constain error according to specified errorMax
    errorIntegral_ += error * deltaT;
    errorIntegral_ = max(min(errorIntegral_, integralErrorMax_), -integralErrorMax_);
    const scalar errorDifferential = (error - oldError_);
    oldError_ = error;

    // Calculate increased output RPS value

    scalar increasedOutputSignal = P_*error + I_*errorIntegral_ + D_*errorDifferential;
    const scalar deltaMax = abs(deltaT*cRate_);
    if(abs(increasedOutputSignal)>= deltaMax)
    {
      increasedOutputSignal = increasedOutputSignal/(abs(increasedOutputSignal)+VSMALL)*deltaMax;
    }
    //outputSignal_ += increasedOutputSignal;
    outputSignal_ = increasedOutputSignal/deltaT;
    // Return result within defined regulator saturation: outputMax_ and outputMin_
    return max(min(outputSignal_, outputMax_), outputMin_);
}

void coursekeepingControl::write(Ostream &os) const
{
    controlMethod::write(os);
    os.beginBlock("parameters");
    os.writeEntry("Kp", P_);
    os.writeEntry("Ti", I_);
    os.writeEntry("Td", D_);
    os.writeEntryIfDifferent("outputMax", 1., outputMax_);
    os.writeEntryIfDifferent("outputMin", 0., outputMin_);
    os.writeEntryIfDifferent("errMax", VGREAT, errorMax_);
    os.writeEntryIfDifferent("errIntegMax", VGREAT, integralErrorMax_);
    os.endBlock();
}

void coursekeepingControl::write(dictionary& dict) const
{
    controlMethod::write(dict);
    dict.add("controllerP", P_);
    dict.add("controllerI", I_);
    dict.add("controllerD", D_);
    dict.add("controllerTarget", cTarget_);
    dict.add("controllerRate", cRate_);
    dict.add("controllerMax", outputMax_);
    dict.add("controllerMin", outputMin_);
    dict.add("oldError", oldError_);
    dict.add("errorIntegral", errorIntegral_);
    dict.add("outputSignal", outputSignal_);
}

