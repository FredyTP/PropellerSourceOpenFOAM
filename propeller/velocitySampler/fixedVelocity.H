#ifndef PROPELLER_VELOCITY_SAMPLER_FIXED_VELOCITY
#define PROPELLER_VELOCITY_SAMPLER_FIXED_VELOCITY
#include "velocitySampler.H"
#include "rotorDiscrete.H"
#include "rotorMesh.H"
#include "dictionary.H"
namespace Foam
{
class fixedVelocity : public velocitySampler
{
public:
    TypeName("fixedVelocity");

    fixedVelocity(const dictionary& dict,const rotorDiscrete* rDiscrete_,const rotorMesh* rMesh_);

    bool read(const dictionary& dict);

    virtual const vectorField& sampleVelocity(const volVectorField& U);

private:

    vector velocity;
};
}

#endif /* PROPELLER_VELOCITY_SAMPLER_FIXED_VELOCITY */


