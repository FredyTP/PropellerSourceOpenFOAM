#ifndef PROPELLER_VELOCITY_SAMPLER_OFFSET_SAMPLER
#define PROPELLER_VELOCITY_SAMPLER_OFFSET_SAMPLER

#include "velocitySampler.H"
#include "rotorDiscrete.H"
#include "rotorMesh.H"
#include "dictionary.H"

namespace Foam
{
class offsetSampler : public velocitySampler
{
public:
    TypeName("offsetSampler");

    offsetSampler(const dictionary& dict,const rotorDiscrete* rDiscrete_,const rotorMesh* rMesh_);

    bool read(const dictionary& dict);

    vector sampleVelocityAt(const volVectorField& U,label i) const;

    bool build();
private:
    scalar offset;
    bool atCellCenter;
    List<label> cellToSample;
    List<point> posToSample;
};

}

#endif /* PROPELLER_VELOCITY_SAMPLER_OFFSET_SAMPLER */


