#ifndef PROPELLER_VELOCITY_SAMPLER_DOMAIN_SAMPLER
#define PROPELLER_VELOCITY_SAMPLER_DOMAIN_SAMPLER

#include "diskSampler.H"
#include "rotorGrid.H"
#include "rotorFvMeshSel.H"
#include "dictionary.H"
#include "cellPointWeight.H"
#include "vector.H"

namespace Foam
{

/**
 * @brief Sampled velocity class from fv vector field solution
 *  
 * This class samples fvSolution to specified points, if no offset is provided,
 * the velocity is sampled at rotor cells.
 * The offset is from each point, with the direction of the rotor normal (upstream)
 * 
 * If velocity is sampled at cell center, its cheaper but less exact, and the value 
 * is constant within a fvCell.
 * 
 * If velocity is not sampled at cell center, the value is interpolated with a tetraedron
 * discrete algorithm
 * @example
 * 
 * Example of dicts provided to this class
 * 
 *   velocitySampler
 *   {
 *   
 *       type            domainSampler;
 *       offset          0.2;
 *       atCellCenter    true;
 *   }
 */
template<class fType>
class domainSampler : public diskSampler<fType>
{
public:
    TypeName("domainSampler");

    domainSampler(const dictionary& dict,const rotorGrid* rGrid,const rotorFvMeshSel* rMesh);

    bool read(const dictionary& dict);

    virtual const Field<fType>& sampleField(const GeometricField<fType, fvPatchField, volMesh>& U);

    void build();

    void writeSampled(const word& name);

    bool isDirectSample();
private:
    scalar offset;
    scalar scale;
   
    List<label> cellToSample;
    List<point> posToSample;
    List<autoPtr<cellPointWeight>> cellWeights;
};

}

#include "domainSampler.C"

#endif /* PROPELLER_VELOCITY_SAMPLER_DOMAIN_SAMPLER */


