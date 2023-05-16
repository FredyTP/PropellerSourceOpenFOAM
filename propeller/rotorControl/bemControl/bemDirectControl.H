#ifndef PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_DIRECT_CONTROL
#define PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_DIRECT_CONTROL

namespace Foam
{

class bemDirectControl : public bemControl
{
public:
    TypeName("bemDirectControl");

    bemDirectControl(const dictionary& dict);
};


    scalar collectivePitch, sinPitch, cosPitch;
    scalar sinAzimuth, cosAzimuth;
    scalar flapping, sinFlapping, cosFlapping;

}

#endif /* PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_DIRECT_CONTROL */
