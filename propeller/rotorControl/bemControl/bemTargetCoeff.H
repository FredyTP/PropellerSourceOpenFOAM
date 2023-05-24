#ifndef PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_TARGET_COEFF
#define PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_TARGET_COEFF

#include "bemControl.H"


namespace Foam
{

class bemTargetCoeff : public bemControl
{
public:
    TypeName("bemTargetCoeff");

    bemTargetCoeff(const dictionary& dict, const bladeElementModel& bem);

    void read(const dictionary& dict);

    void correctControl(const vectorField& U, const scalarField* rhoField);

    scalar getAzimuth(scalar azimuth0) const;
    scalar getPitch(scalar azimuth) const;
    scalar getFlapping(scalar azimuth) const;

    scalar getAzimuthDot(scalar azimuth) const;
    scalar getPitchDot(scalar azimuth) const;
    scalar getFlappingDot(scalar azimuth) const;

    vector calcForces(const vectorField &U, const scalarField *rhoField) const;
    virtual scalar getAngularVelocity() const;

private:


    const bladeElementModel& bem_;

    //- Number of iterations between calls to 'correct'
    label calcFrequency_;

    //- Target coefficient vector (thrust force, roll moment, pitch moment)
    vector target_;

    //- Pitch angles (collective, roll, pitch) [rad]
    vector theta_;

    //- Maximum number of iterations in trim routine
    label nIter_;

    //- Convergence tolerance
    scalar tol_;

    //- Under-relaxation coefficient
    scalar relax_;

    //- Perturbation angle used to determine jacobian
    scalar dTheta_;

    scalar angularVelocity_;
};




}


#endif /* PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_TARGET_COEFF */
