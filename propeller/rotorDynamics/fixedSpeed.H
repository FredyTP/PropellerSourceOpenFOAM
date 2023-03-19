#ifndef PROPELLER_ROTOR_DYNAMICS_FIXED_SPEED
#define PROPELLER_ROTOR_DYNAMICS_FIXED_SPEED

#include "rotorDynamics.H"
#include "dictionary.H"


namespace Foam
{

class fixedSpeed : public rotorDynamics
{
public:

    TypeName("fixedSpeed");

    fixedSpeed(const dictionary& dict);

    virtual void integrate(scalar aeroMoment,scalar dt)
    {
        Info<<"time step: "<<dt<<endl;
        torque_  = (-aeroMoment + omega_ * viscousDisipation_);
    }
private:


};





}


















#endif //PROPELLER_ROTOR_DYNAMICS_FIXED_SPEED