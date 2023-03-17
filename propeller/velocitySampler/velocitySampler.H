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
     * Return sampled velocities on rotor ordered cells positions(rotor cells goes from 0 to N)
     * Do not get confused by celli which usually means fvmesh cell index
    */
    virtual const vectorField& sampleVelocity(const volVectorField& U) = 0;

    virtual void writeSampled(const word& name);
    
    static autoPtr<velocitySampler> New(const dictionary& dict,const rotorDiscrete* rDiscrete_,const rotorMesh* rMesh_);

protected:
    velocitySampler(const rotorDiscrete* rDiscrete_,const rotorMesh* rMesh_);
    const rotorMesh* rMesh;
    const rotorDiscrete* rDiscrete;
    vectorField sampledVel;
};


}

#endif /* PROPELLER_VELOCITY_SAMPLER_VELOCITY_SAMPLER */


