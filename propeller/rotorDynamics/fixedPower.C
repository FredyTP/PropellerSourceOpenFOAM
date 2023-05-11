#include "fixedPower.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    //set and define the type name "typeName"
    defineTypeNameAndDebug(fixedPower,0);

    //Add to run time table to dynamically select the class based on
    //the dictionary airfoilModel atribute

    addToRunTimeSelectionTable(rotorDynamics,fixedPower,dictionary);

    fixedPower::fixedPower()
    {
    }

    fixedPower::fixedPower(const dictionary &dict)
    {
        bool isRpm = false;
        scalar rpm;
        isRpm = dict.readIfPresent("rpm",rpm);

        if(isRpm)
        {
            omega_ = rpmToRad_s(rpm);
        }
        else
        {
            omega_ = dict.getOrDefault<scalar>("angularRate",0.0);
        }

        inertia_ = dict.getOrDefault<scalar>("angularRate",0.0);
        torque_ = dict.getOrDefault<scalar>("torque",0.0);

        //TODO: finish this ...

    }


}
