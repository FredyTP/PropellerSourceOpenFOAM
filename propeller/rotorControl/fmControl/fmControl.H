#ifndef PROPELLER_ROTOR_CONTROL_FM_CONTROL_BEM_CONTROL
#define PROPELLER_ROTOR_CONTROL_FM_CONTROL_BEM_CONTROL
#ifndef PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_CONTROL
#define PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_CONTROL

#include "rotorControl.H"
#include "unitConversion.H"
#include "vectorField.H"
#include "scalarField.H"

namespace Foam
{

class forceModel;
class fmControl : public rotorControl
{
public:
    TypeName("fmControl");
    
    declareRunTimeSelectionTable
    (
        autoPtr,
        fmControl,
        dictionary,
        (
            const dictionary& dict,
            const forceModel& model
        ),
        (dict,model)
    );

    fmControl(const dictionary& dict, const forceModel& fmModel);

    scalar readAngularVelocity(const dictionary& dict);

    virtual void correctControl(const vectorField& U, const scalarField* rhoField){}

    virtual scalar getAngularVelocity() {return 10;};

    scalar getJ();


    static autoPtr<fmControl> New(const dictionary& dict, const forceModel& bem);
private:
    scalar J=0.5;
};

}


#endif /* PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_CONTROL */


#endif /* PROPELLER_ROTOR_CONTROL_FM_CONTROL_BEM_CONTROL */
