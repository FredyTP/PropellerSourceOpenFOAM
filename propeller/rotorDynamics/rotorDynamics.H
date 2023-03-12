#ifndef PROPELLER_ROTOR_DYNAMICS_ROTOR_DYNAMICS
#define PROPELLER_ROTOR_DYNAMICS_ROTOR_DYNAMICS

#include "scalar.H"

namespace Foam
{

class rotorDynamics
{
public:
    rotorDynamics()
    {

    }

    /**
     * Integrates for delta t to obtain new angular velocity
     * Positive moment is clockwise
     * Simple order 1 algorithm used
    */
    virtual void integrate(scalar aeroMoment,scalar dt)
    {
        scalar dwdt = (aeroMoment + par)/Izz;
        omega += dwdt*dt;
    }

    scalar angularVelocity() const
    {
        return omega;
    }
private:
    scalar omega=10; // rad/s
    scalar par=200; // Nm
    scalar Izz = 1;// kgm
};

}

#endif /* PROPELLER_ROTOR_DYNAMICS_ROTOR_DYNAMICS */
