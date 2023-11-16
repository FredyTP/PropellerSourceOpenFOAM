#ifndef PROPELLER_VELOCITY_SAMPLER_FILE_SAMPLER
#define PROPELLER_VELOCITY_SAMPLER_FILE_SAMPLER


#include "domainSampler.H"
#include "rotorGrid.H"
#include "rotorFvMeshSel.H"
#include "dictionary.H"
namespace Foam
{

/**
 * @brief Fixed velocity class, independent of fv Velocity field
 * 
 * Same as offsetSampler but the sampled velocity is from a specified file
 * 
 * @example
 * 
 *   velocitySampler
 *   {
 *   
 *       type            fileSampler;
 *       file            somefile;
 *       offset          0.2;
 *       atCellCenter    true;
 *   }
 * 
 * 
 */

template<class fType>
class fileSampler : public domainSampler<fType>
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
    fileSampler(const dictionary& dict,const rotorGrid* rGrid,const rotorFvMeshSel* rMesh);

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
    virtual const Field<fType>& sampleField(const GeometricField<fType, fvPatchField, volMesh>& U);
private:
    autoPtr<GeometricField<fType, fvPatchField, volMesh>> uFromFile_;
};

}

#include "fileSampler.C"

#endif /* PROPELLER_VELOCITY_SAMPLER_FILE_SAMPLER */



