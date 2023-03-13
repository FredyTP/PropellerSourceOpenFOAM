#ifndef PROPELLER_ROTOR_DYNAMICS_ROTOR_DYNAMICS
#define PROPELLER_ROTOR_DYNAMICS_ROTOR_DYNAMICS

#include "scalar.H"
#include "dictionary.H"
#include "mathematicalConstants.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class rotorDynamics
{
public:
    TypeName("rotorDynamics");

    declareRunTimeSelectionTable
    (
        autoPtr,
        rotorDynamics,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    rotorDynamics();

    virtual void integrate(scalar aeroMoment,scalar dt) = 0;


    static Foam::autoPtr<Foam::rotorDynamics> New(const dictionary& dict);

    scalar angularVelocity() const
    {
        return omega_;
    }
    scalar rpm() const
    {
        return omega_ * 30.0/Foam::constant::mathematical::pi;
    }
    /**
     * Engine torque
    */
    scalar torque() const
    {
        return torque_;
    }

    /**
     * 
    */
    scalar power() const
    {
        return torque_*omega_;
    }

protected:
    scalar omega_; // rad/s
    scalar torque_; // Nm
    scalar inertia_;// kgm
    scalar viscousDisipation_; //kgm/s
};

}

#endif /* PROPELLER_ROTOR_DYNAMICS_ROTOR_DYNAMICS */
