#ifndef PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_DIRECT_CONTROL
#define PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_DIRECT_CONTROL

#include "bemControl.H"

namespace Foam
{

class bemDirectControl : public bemControl
{
public:
    TypeName("bemDirectControl");

    bemDirectControl(const dictionary& dict);

    scalar getAzimuth(scalar azimuth0) const;
    scalar getPitch(scalar azimuth) const;
    scalar getFlapping(scalar azimuth) const;

    scalar getAzimuthDot(scalar azimuth, scalar angularVelocity);
    scalar getPitchDot(scalar azimuth, scalar angularVelocity);
    scalar getFlappingDot(scalar azimuth, scalar angularVelocity);

    vector getBladeLocalVel(scalar azimuth, scalar angularVelocity, scalar radius);

private:
    scalar collectivePitch, sinPitch, cosPitch;
    scalar sinAzimuth, cosAzimuth;
    scalar flapping, sinFlapping, cosFlapping;
};




}

#endif /* PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_DIRECT_CONTROL */
