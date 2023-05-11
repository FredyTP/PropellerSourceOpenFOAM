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

    virtual ~rotorDynamics() = default;

    virtual void integrate(scalar aeroMoment,scalar dt) = 0;


    static Foam::autoPtr<Foam::rotorDynamics> New(const dictionary& dict);

    scalar theta()
    {
        return theta_;
    }
    scalar angularVelocity() const
    {
        return omega_;
    }
    scalar rpm() const
    {
        return rad_sToRpm(omega_);
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

    inline static scalar rpmToRad_s(scalar rpm){ return rpm * rpmToRadFactor;}
    inline static scalar rad_sToRpm(scalar rad_s){ return rad_s * radToRpmFactor;}
protected:
    scalar omega_; // rad/s
    scalar torque_; // Nm
    scalar inertia_;// kgm
    scalar viscousDisipation_; //kgm/s
    scalar theta_=0;

    static constexpr scalar rpmToRadFactor = Foam::constant::mathematical::pi/30.0;
    static constexpr scalar radToRpmFactor = 30.0/Foam::constant::mathematical::pi;
};

}

#endif /* PROPELLER_ROTOR_DYNAMICS_ROTOR_DYNAMICS */


