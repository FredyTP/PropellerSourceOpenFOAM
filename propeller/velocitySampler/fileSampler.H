#ifndef PROPELLER_VELOCITY_SAMPLER_FILE_SAMPLER
#define PROPELLER_VELOCITY_SAMPLER_FILE_SAMPLER

#include "velocitySampler.H"
#include "domainSampler.H"
#include "rotorDiscrete.H"
#include "rotorFvMeshSel.H"
#include "dictionary.H"
namespace Foam
{

/**
 * @brief Fixed velocity class, independent of fv Velocity field
 * 
 * Same as offsetSampler but the sampled velocity is from a specified file
 * 
 * Example of dicts provided to this class
 * 
 *   velocitySampler
 *   {
 *   
 *       type            fileSampler;
 *       offset          0.2;
 *       atCellCenter    true;
 *   }
 * 
 * 
 */
class fileSampler : public domainSampler
{
public:
    TypeName("fileSampler");

    /**
     * @brief Construct a new fixed Velocity object
     * 
     * @param dict 
     * @param rDiscrete_ 
     * @param rMesh_ 
     */
    fileSampler(const dictionary& dict,const rotorDiscrete* rDiscrete_,const rotorFvMeshSel* rMesh_);

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
    autoPtr<volVectorField> uFromFile_;
    vector velocity;
};
}

#endif /* PROPELLER_VELOCITY_SAMPLER_FILE_SAMPLER */



