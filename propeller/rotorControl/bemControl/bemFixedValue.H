#ifndef PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_FIXED_VALUE
#define PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_FIXED_VALUE

#include "bemControl.H"

namespace Foam
{

class bemFixedValue : public bemControl
{
public:
    TypeName("bemFixedValue");

    bemFixedValue(const dictionary& dict, const bladeElementModel& bem);

    void correctControl(const vectorField& U, const scalarField* rhoField);

    scalar getAzimuth(scalar azimuth0) const;
    scalar getPitch(scalar azimuth) const;
    scalar getFlapping(scalar azimuth) const;

    scalar getAzimuthDot(scalar azimuth) const;
    scalar getPitchDot(scalar azimuth) const;
    scalar getFlappingDot(scalar azimuth) const;

    scalar getAngularVelocity() const;
private:
    scalar collectivePitch, sinPitch, cosPitch;
    scalar sinAzimuth, cosAzimuth;
    scalar flapping, sinFlapping, cosFlapping;
    
    scalar angularVelocity_;
};




}

#endif /* PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_FIXED_VALUE */
