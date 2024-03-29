#include "rotorControl.H"
#include "runTimeSelectionTables.H"
#include "unitConversion.H"

namespace Foam
{

scalar rotorControl::readAngularVelocity(const dictionary &dict)
{
    scalar rpm;

    if(dict.readIfPresent("rpm",rpm))
    {
        return rpmToRads(rpm);
    }
    else
    {
        return dict.get<scalar>("angularVelocity");
    }
}

}