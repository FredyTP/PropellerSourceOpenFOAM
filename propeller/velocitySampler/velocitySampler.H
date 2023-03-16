#ifndef PROPELLER_VELOCITY_SAMPLER_VELOCITY_SAMPLER
#define PROPELLER_VELOCITY_SAMPLER_VELOCITY_SAMPLER
#include "dictionary.H"
#include "autoPtr.H"
#include "rotorDiscrete.H"
#include "rotorMesh.H"

namespace Foam
{
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
            const rotorMesh* rMesh_
        ),
        (dict,rDiscrete_,rMesh_)
    );

    /**
     * Return sampled velocity by rotor ordered cell at index i (rotor cells goes from 0 to N)
     * Do not get confused by celli which usually means fvmesh cell index
    */
    virtual vector sampleVelocityAt(const volVectorField& U,label i) const = 0;


    
    static autoPtr<velocitySampler> New(const dictionary& dict,const rotorDiscrete* rDiscrete_,const rotorMesh* rMesh_);

protected:
    velocitySampler(const rotorDiscrete* rDiscrete_,const rotorMesh* rMesh_);
    const rotorMesh* rMesh;
    const rotorDiscrete* rDiscrete;
};


}

#endif /* PROPELLER_VELOCITY_SAMPLER_VELOCITY_SAMPLER */


