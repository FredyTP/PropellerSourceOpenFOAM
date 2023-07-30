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

/**
 * @brief Base abstract class of rotor control with some helper methods and pure virtual method to implement
 * 
 */
class rotorControl
{
public:
    static scalar readAngularVelocity(const dictionary& dict);
    virtual void correctControl(const vectorField& U, const scalarField* rhoField) = 0;
    virtual ~rotorControl() = default;
};

/**
 * @brief Templated base classs including acount for output and control vars via typedef of template argument
 * 
 * @tparam model used to control with the variable definitions
 */
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
