#ifndef PROPELLER_VELOCITY_SAMPLER_AZIMUTAL_MEAN
#define PROPELLER_VELOCITY_SAMPLER_AZIMUTAL_MEAN

#include "domainSampler.H"
#include "rotorGrid.H"
#include "polarGrid.H"
#include "rotorFvMeshSel.H"
#include "dictionary.H"
namespace Foam
{

/**
 * @brief This class extends domainSampler making an azimutal mean of the cells
 * This class only works with polarGrid
 * @example
 * 
 * Example of dicts provided to this class
 * 
 *   velocitySampler
 *   {
 *   
 *       type            azimutalMean;
 *       offset          0.2;
 *       scale 1;
 *       atCellCenter    true;
 *   }
 */
class azimutalMean : public domainSampler
{
public:
    TypeName("azimutalMean");

    azimutalMean(const dictionary& dict,const rotorGrid* rGrid,const rotorFvMeshSel* rMesh);

    virtual const vectorField& sampleVelocity(const volVectorField& U);

private:

    const polarGrid* polarGrid_;

};

}

#endif /* PROPELLER_VELOCITY_SAMPLER_AZIMUTAL_MEAN */


