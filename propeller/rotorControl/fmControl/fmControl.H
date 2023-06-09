#ifndef PROPELLER_ROTOR_CONTROL_FM_CONTROL_FM_CONTROL
#define PROPELLER_ROTOR_CONTROL_FM_CONTROL_FM_CONTROL

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

    virtual void correctControl(const vectorField& U, const scalarField* rhoField) = 0;

    virtual scalar getAngularVelocity() const = 0;

    virtual scalar getJ() const = 0;

    static autoPtr<fmControl> New(const dictionary& dict, const forceModel& bem);

    scalar getJFromOmega(scalar omega, scalar speed);
    scalar getOmegaFromJ(scalar J, scalar speed);
    
protected:
    
    const forceModel& forceModel_;
};

}


#endif /* PROPELLER_ROTOR_CONTROL_FM_CONTROL_FM_CONTROL */
