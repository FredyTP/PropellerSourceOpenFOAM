#ifndef PROPELLER_ROTOR_CONTROL_FM_CONTROL_FM_TARGET_VALUE
#define PROPELLER_ROTOR_CONTROL_FM_CONTROL_FM_TARGET_VALUE

#include "fmControl.H"
#include "unitConversion.H"
#include "vectorField.H"
#include "scalarField.H"

namespace Foam
{

class forceModel;
class fmTargetValue : public fmControl
{
public:

    TypeName("fmTargetValue");

    fmTargetValue(const dictionary& dict, const forceModel& fmModel);

    scalar readAngularVelocity(const dictionary& dict);

    virtual void correctControl(const vectorField& U, const scalarField* rhoField);

    virtual scalar getAngularVelocity() const;

    virtual scalar getJ() const;

 
private:
    std::function<scalarField(scalarField)> forceFunction(const vectorField& U, const scalarField* rhoField);
    std::function<scalarField(scalarField)> torqueFunction(const vectorField& U, const scalarField* rhoField);

    bool targetThrust_;
    scalar thrust_;
    scalar torque_;

    scalar J_;
    scalar omega_;

};

}




#endif /* PROPELLER_ROTOR_CONTROL_FM_CONTROL_FM_TARGET_VALUE */
