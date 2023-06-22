#ifndef PROPELLER_ROTOR_CONTROL_ROTOR_CONTROL
#define PROPELLER_ROTOR_CONTROL_ROTOR_CONTROL
#include "dictionary.H"
#include "autoPtr.H"
#include "scalar.H"
#include "vector.H"
#include <type_traits>
#include "propellerModel.H"
#include "EnumField.H"
namespace Foam
{

class rotorControl
{
public:
    static scalar readAngularVelocity(const dictionary& dict);
    virtual ~rotorControl() = default;
};

template<class model>
class rotorControlTemplate : public rotorControl
{
public:

    //static_assert(std::is_base_of<propellerModel,model>::value,"rotorControl must be templated with a propellerModel");

    typedef typename model::outputVar outputVar;
    typedef typename model::controlVar controlVar;

    typedef util::EnumField<outputVar,model::N_output> outputVarType;
    typedef util::EnumField<controlVar,model::N_control> controlVarType;

    
    
};

}

#endif /* PROPELLER_ROTOR_CONTROL_ROTOR_CONTROL */
