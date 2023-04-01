#ifndef PROPELLER_VELOCITY_SAMPLER_FIXED_VELOCITY
#define PROPELLER_VELOCITY_SAMPLER_FIXED_VELOCITY
#include "velocitySampler.H"
#include "rotorDiscrete.H"
#include "rotorFvMeshSel.H"
#include "dictionary.H"
namespace Foam
{

/**
 * @brief Fixed velocity class, independent of fv Velocity field
 * 
 * Velocity specification can be surface normal, or absolute velocity.
 * If velocity is surface normal, the velocity parameter is an scalar.
 * -Positive velocity values are towards the rotor, oposite to rotor normal vector.
 * 
 * If velocity is absolute, the velocity parameter is a vector.
 * @example
 * 
 * Example of dicts provided to this class
 * 
 *   velocitySampler
 *   {
 *   
 *       type            fixedVelocity;
 *       normal          true;
 *       velocity        10;
 *   }
 * 
 *   velocitySampler
 *   {
 *   
 *       type            fixedVelocity;
 *       normal          false;
 *       velocity        (1, -10 , 15);
 *   }
 * 
 * 
 */
class fixedVelocity : public velocitySampler
{
public:
    TypeName("fixedVelocity");

    /**
     * @brief Construct a new fixed Velocity object
     * 
     * @param dict 
     * @param rDiscrete_ 
     * @param rMesh_ 
     */
    fixedVelocity(const dictionary& dict,const rotorDiscrete* rDiscrete_,const rotorFvMeshSel* rMesh_);

    /**
     * @brief Read object from dictionary
     * 
     * @param dict 
     * @return true 
     * @return false 
     */
    bool read(const dictionary& dict);

    /**
     * @brief Just return the constant velocity field stored in base class with no additional computations
     * 
     * @param U 
     * @return const vectorField& 
     */
    virtual const vectorField& sampleVelocity(const volVectorField& U);

private:
    //- Value of the constant velocity field
    vector velocity;
};
}

#endif /* PROPELLER_VELOCITY_SAMPLER_FIXED_VELOCITY */


