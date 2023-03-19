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
        

    }


}

