#ifndef PROPELLER_VELOCITY_SAMPLER_VELOCITY_SAMPLER
#define PROPELLER_VELOCITY_SAMPLER_VELOCITY_SAMPLER
#include "dictionary.H"
#include "autoPtr.H"
#include "rotorDiscrete.H"
#include "rotorFvMeshSel.H"

namespace Foam
{

/**
 * @brief Velocity sample base class, used to provide velocity data to
 * propeller model
 * 
 * 
 * @example
 * 
 * Example of dicts provided to this class
 * 
 * -Fixed velocity sampler
 * 
 *   velocitySampler
 *   {
 *   
 *       type            fixedVelocity;
 *       normal          true;
 *       velocity        10;
 *   }
 */
class velocitySampler
{
public:
    TypeName("velocitySampler");
    declareRunTimeSelectionTable
    (
        autoPtr,
        velocitySampler,
        dictionary,
        (
            const dictionary& dict,
            const rotorDiscrete* rDiscrete_,
            const rotorFvMeshSel* rMesh_
        ),
        (dict,rDiscrete_,rMesh_)
    );

    /**
     * @brief Return sampled velocities on rotor discrete positions(rotor positions goes from 0 to N)
     * Do not get confused by celli which usually means fvmesh cell index
     * 
     * @param U fv velocity field
     * @return const vectorField& 
     */
    virtual const vectorField& sampleVelocity(const volVectorField& U) = 0;

    /**
     * @brief Write field of sampled cells used to interpolate velocity field
     * 
     * @param name Propeller's name for output, adding ":sampledCells" to that name.
     *  Ex: propeller1:sampledCells
     */
    virtual void writeSampled(const word& name);
    
    /**
     * @brief Builds a velocitySampler from its dictionary, rotorDiscrete data and fv Mesh selection
     * 
     * @param dict velocity sampler dict
     * @param rDiscrete_ rotor discretization
     * @param rMesh_ fvMesh selection for rotor
     * @return autoPtr<velocitySampler>: safe pointer to created sampler
     */
    static autoPtr<velocitySampler> New(const dictionary& dict,const rotorDiscrete* rDiscrete_,const rotorFvMeshSel* rMesh_);

protected:

    /**
     * @brief Construct a new velocity Sampler object, only available to build from derived clases
     * 
     * @param rDiscrete_ 
     * @param rMesh_ 
     */
    velocitySampler(const rotorDiscrete* rDiscrete_,const rotorFvMeshSel* rMesh_);

    const rotorFvMeshSel* rMesh;
    const rotorDiscrete* rDiscrete;

    /**
     * @brief velocity sampled vector, to return by reference
     * or be more efficient when constant velocity is required
     * 
     */
    vectorField sampledVel;
};


}

#endif /* PROPELLER_VELOCITY_SAMPLER_VELOCITY_SAMPLER */


