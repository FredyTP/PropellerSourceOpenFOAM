#ifndef PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_FIXED_CONTROL
#define PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_FIXED_CONTROL

#include "bemControl.H"

namespace Foam
{


/**
 * @brief Class for fixedc control values of blade element method via dictionary
 * 
 */
class bemFixedControl : public bemControl
{
public:
    TypeName("fixedControl");

    bemFixedControl(const dictionary& dict, const bladeElementModel& bem);

    void correctControl(const vectorField& U, const scalarField* rhoField);

    scalar getAzimuth(scalar azimuth0) const;
    scalar getPitch(scalar azimuth) const;
    scalar getFlapping(scalar azimuth) const;

    scalar getAzimuthDot(scalar azimuth) const;
    scalar getPitchDot(scalar azimuth) const;
    scalar getFlappingDot(scalar azimuth) const;

    scalar getAngularVelocity() const;
private:
    //- Control var (angularVelocity collective roll pitch) [rad/s rad rad rad]
    controlVarType control_;

    //- Lagging parameters
    scalar sinAzimuth, cosAzimuth;

    //- Flapping parameters
    scalar flapping, sinFlapping, cosFlapping;
};




}

#endif /* PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_FIXED_CONTROL */
