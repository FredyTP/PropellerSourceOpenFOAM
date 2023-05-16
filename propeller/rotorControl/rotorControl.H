#ifndef PROPELLER_ROTOR_CONTROL_ROTOR_CONTROL
#define PROPELLER_ROTOR_CONTROL_ROTOR_CONTROL
#include "dictionary.H"
#include "autoPtr.H"

namespace Foam
{
class rotorControl
{
public:
    TypeName("rotorControl");
    declareRunTimeSelectionTable
    (
        autoPtr,
        rotorControl,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );
};
}

#endif /* PROPELLER_ROTOR_CONTROL_ROTOR_CONTROL */
