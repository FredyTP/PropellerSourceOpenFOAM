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

    fmControl(const dictionary& dict);

    scalar readAngularVelocity(const dictionary& dict);

    virtual void correctControl(const vectorField& U, const scalarField* rhoField) = 0;

    virtual scalar getAzimuth(scalar azimuth0) const = 0;
    virtual scalar getPitch(scalar azimuth) const = 0;
    virtual scalar getFlapping(scalar azimuth) const = 0;

    virtual scalar getAzimuthDot(scalar azimuth) const = 0;
    virtual scalar getPitchDot(scalar azimuth) const = 0;
    virtual scalar getFlappingDot(scalar azimuth) const = 0;

    virtual scalar getAngularVelocity() const = 0;

    //virtual scalar getAngularVelocity() = 0;
    vector getBladeLocalVel(scalar azimuth, scalar radius) const;


    static autoPtr<fmControl> New(const dictionary& dict, const forceModel& bem);
};

}


#endif /* PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_CONTROL */


#endif /* PROPELLER_ROTOR_CONTROL_FM_CONTROL_BEM_CONTROL */
