#ifndef PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_CONTROL
#define PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_CONTROL

#include "rotorControl.H"

namespace Foam
{
class bemControl : public rotorControl
{
public:
    TypeName("bemControl");
    
    bemControl(const dictionary& dict);

};
}


#endif /* PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_CONTROL */
