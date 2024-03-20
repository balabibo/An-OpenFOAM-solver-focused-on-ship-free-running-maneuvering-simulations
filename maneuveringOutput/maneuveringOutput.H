/*---------------------------------------------------------------------------*\
Class
    Foam::maneuveringOutput

Description
    Test

\*---------------------------------------------------------------------------*/

#ifndef Foam_maneuveringOutput_H
#define Foam_maneuveringOutput_H

#include "fvCFD.H"
#include "maneuveringInput.H"
#include "controlMethod.H"

//- Class defining different sorts of maneuvering motions
class maneuveringOutput
{
public:
    maneuveringOutput() = delete;
    ~maneuveringOutput() = default;

    // Initialize from mesh and dynamicMeshDict dictionary
    maneuveringOutput(const polyMesh &mesh, const dictionary &dict);

    //- Initialize from mesh, dynamicMeshDict dictionary and control state dictionary
    //maneuveringOutput(const polyMesh &mesh, const dictionary &dict, const dictionary& stateDict);

    //- output the results after maneuvering controller
    scalar output(const vector2D& input);
     
    
    //void write(Ostream& os, const word dictName = "maneuveringOutput") const;
    //- Write to to stream
    void write(Ostream& os) const;

    //- Write to dictionary
    void write(dictionary& dict) const;

    //- reture shared_ptr pointing Input
    const std::shared_ptr<maneuveringInput> mInput();

private:

    //- To acquire time
    const polyMesh &mesh_;

    //- List of references of maneuvering controller
    const std::shared_ptr<maneuveringInput> maneuveringInput_;
    
    //- List of control algorithms
    std::shared_ptr<controlMethod> controlMethod_;

    //- Time index of present update
    label timeIndex_;
};



#endif

