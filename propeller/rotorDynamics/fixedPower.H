#ifndef PROPELLER_ROTOR_DYNAMICS_FIXED_POWER
#define PROPELLER_ROTOR_DYNAMICS_FIXED_POWER

#include "scalar.H"
#include "rotorDynamics.H"

namespace Foam
{

class fixedPower : public rotorDynamics
{
public:
    TypeName("fixedPower");
    
    fixedPower();

    fixedPower(const dictionary& dict);
    /**
     * Integrates for delta t to obtain new angular velocity
     * Positive moment is clockwise
     * Simple order 1 algorithm used
    */
    virtual void integrate(scalar aeroMoment,scalar dt)
    {
        if(inertia_ !=0)
        {
            scalar dwdt = (aeroMoment + torque_ - viscousDisipation_*omega_)/inertia_;
            omega_ += dwdt*dt;
        }
        else
        {
            //Solve equation using Ma(w)
            omega_ = 0;
        }
    }

private:
    
};

}

#endif /* PROPELLER_ROTOR_DYNAMICS_ROTOR_DYNAMICS */

