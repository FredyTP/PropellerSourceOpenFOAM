#ifndef PROPELLER_VELOCITY_SAMPLER_FIXED_VALUE
#define PROPELLER_VELOCITY_SAMPLER_FIXED_VALUE
#include "diskSampler.H"
#include "rotorGrid.H"
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
 *   diskSampler
 *   {
 *   
 *       type            fixedValue;
 *       normal          true;
 *       velocity        10;
 *   }
 * 
 *   diskSampler
 *   {
 *   
 *       type            fixedValue;
 *       normal          false;
 *       velocity        (1, -10 , 15);
 *   }
 * 
 * 
 */
template<class fType>
class fixedValue : public diskSampler<fType>
{
public:
    TypeName("fixedValue");

    /**
     * @brief Construct a new fixed Velocity object
     * 
     * @param dict 
     * @param rDiscrete_ 
     * @param rMesh_ 
     */
    fixedValue(const dictionary& dict,const rotorGrid* rGrid,const rotorFvMeshSel* rMesh);

    /**
     * @brief Read object from dictionary
     * 
     * @param dict 
     * @return true 
     * @return false 
     */
    bool read(const dictionary& dict);

    /**
     * @brief Just return the constant field stored in base class with no additional computations
     * 
     * @param U 
     * @return const Field<fType>& 
     */
    virtual const Field<fType>& sampleField(const GeometricField<fType, fvPatchField, volMesh>& U);

private:
    //- Value of the constant field
    fType fieldValue_;
};
}

#include "fixedValue.C"

#endif /* PROPELLER_VELOCITY_SAMPLER_FIXED_VELOCITY */


