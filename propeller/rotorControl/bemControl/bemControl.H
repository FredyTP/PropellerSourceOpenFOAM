#ifndef PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_CONTROL
#define PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_CONTROL

#include "rotorControl.H"
#include "unitConversion.H"
#include "vectorField.H"
#include "scalarField.H"
#include "bladeElementModel.H"
namespace Foam
{


class bemControl : public rotorControlTemplate<bladeElementModel>
{
public:
    TypeName("bemControl");
    
    declareRunTimeSelectionTable
    (
        autoPtr,
        bemControl,
        dictionary,
        (
            const dictionary& dict,
            const bladeElementModel& model
        ),
        (dict,model)
    );

    bemControl(const dictionary& dict);


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


    static autoPtr<bemControl> New(const dictionary& dict, const bladeElementModel& bem);

    static scalar PitchFunction(scalar azimuth, scalar pitch0, scalar pitch1c, scalar pitch1s);
    static scalar AzimuthFunction(scalar azimuth0, scalar azimuth1c, scalar azimuth1s);
    static scalar FlappingFunction(scalar azimuth, scalar flapping0, scalar flapping1c, scalar flapping1s);

    static scalar PitchDotFunction(scalar omega, scalar azimuth, scalar pitch0, scalar pitch1c, scalar pitch1s);
    static scalar AzimuthDotFunction(scalar omega, scalar azimuth0, scalar azimuth1c, scalar azimuth1s);
    static scalar FlappingDotFunction(scalar omega, scalar azimuth, scalar flapping0, scalar flapping1c, scalar flapping1s);
};

}


#endif /* PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_CONTROL */
