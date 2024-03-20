/*---------------------------------------------------------------------------*\
Class
    Foam::maneuveringInput

Description
    Part of the maneuveringOutput.C， which is used to read sailing velocity or 
yaw angle of prescribed rigid body.

\*---------------------------------------------------------------------------*/

#ifndef Foam_maneuveringInput_H
#define Foam_maneuveringInput_H

#include "fvCFD.H"

class maneuveringInput
{
public:
    virtual ~maneuveringInput() = default;
    maneuveringInput() = delete;

    // Initilize from mesh and dictionary
    maneuveringInput(const dictionary &dict);

    const word inputName() const;
    const word refBody() const;
    const word actBody() const;
    virtual label inputTypeValue() const = 0;
    virtual label controlType() const = 0;
    

    static std::shared_ptr<maneuveringInput> create(const dictionary &dict);

    enum inputType
    {
       sailing,
       turning,
       zigzag,
       coursekeeping
    };

    static const Enum<inputType> inputTypeNames;
    virtual void write(dictionary&) const;

protected:
      
    const word inputName_;
    const word refBody_;
    const word actBody_;   
          
};

class sailingInput : public maneuveringInput
{
public:
    sailingInput() = delete;
    sailingInput(const dictionary &dict);
    label inputTypeValue() const override;      
    label controlType() const override;
    void write(dictionary &) const override;
        
};

class yawInput : public maneuveringInput
{
public:
    yawInput() = delete;
    yawInput(const dictionary &dict);
    label inputTypeValue() const override;
    label controlType() const override;
    void write(dictionary &) const override;

};


#endif
