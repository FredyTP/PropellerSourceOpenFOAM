#include "fixedSpeed.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    //set and define the type name "typeName"
    defineTypeNameAndDebug(fixedSpeed,0);

    //Add to run time table to dynamically select the class based on
    //the dictionary airfoilModel atribute

    addToRunTimeSelectionTable(rotorDynamics,fixedSpeed,dictionary);

    fixedSpeed::fixedSpeed(const dictionary &dict)
    {
        Info.stream().incrIndent();

        bool isRpm = false;
        scalar rpm;
        isRpm = dict.readIfPresent("rpm",rpm);

        if(isRpm)
        {
            omega_ = rpmToRad_s(rpm);
        }
        else
        {
            dict.readEntry("angularRate",omega_);
        }
        
        Info<<indent
            <<"- Angular Rate: "<<omega_<<" (rad/s) / "
            <<this->rpm()<<" (rpm)"<<endl;
        Info.stream().decrIndent();
    }

    void fixedSpeed::integrate(scalar aeroMoment, scalar dt)
    {
        theta_+= dt*this->angularVelocity();
        if(theta_>constant::mathematical::pi)
        {
            theta_ -=constant::mathematical::twoPi;
        }
        torque_  = (-aeroMoment + omega_ * viscousDisipation_);
    }
}

