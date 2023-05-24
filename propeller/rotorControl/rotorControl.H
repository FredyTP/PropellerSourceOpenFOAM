#ifndef PROPELLER_ROTOR_CONTROL_ROTOR_CONTROL
#define PROPELLER_ROTOR_CONTROL_ROTOR_CONTROL
#include "dictionary.H"
#include "autoPtr.H"
#include "scalar.H"
#include "vector.H"
namespace Foam
{
class rotorControl
{
public:
    TypeName("rotorControl");

    
    virtual ~rotorControl() = default;
};

}

#endif /* PROPELLER_ROTOR_CONTROL_ROTOR_CONTROL */
