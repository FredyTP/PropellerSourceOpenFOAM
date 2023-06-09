#ifndef PROPELLER_ROTOR_CONTROL_FM_CONTROL_FM_FIXED_VALUE
#define PROPELLER_ROTOR_CONTROL_FM_CONTROL_FM_FIXED_VALUE


#include "fmControl.H"
#include "unitConversion.H"
#include "vectorField.H"
#include "scalarField.H"

namespace Foam
{

class forceModel;
class fmFixedValue : public fmControl
{
public:
    TypeName("fmFixedValue");

    fmFixedValue(const dictionary& dict, const forceModel& fmModel);

    scalar readAngularVelocity(const dictionary& dict);

    virtual void correctControl(const vectorField& U, const scalarField* rhoField);

    virtual scalar getAngularVelocity() const;

    virtual scalar getJ() const;

private:
    
    bool fixedJ_;
    scalar J_;
    scalar omega_;

};

}




#endif /* PROPELLER_ROTOR_CONTROL_FM_CONTROL_FM_FIXED_VALUE */
