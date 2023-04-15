#include "fixedVelocity.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam
{
 
    defineTypeNameAndDebug(fixedVelocity,0);
    addToRunTimeSelectionTable(velocitySampler,fixedVelocity, dictionary);


fixedVelocity::fixedVelocity(const dictionary& dict,const rotorDiscrete* rDiscrete_,const rotorFvMeshSel* rMesh_)
    : velocitySampler(rDiscrete_,rMesh_)
{
    this->read(dict);
}

bool fixedVelocity::read(const dictionary &dict)
{
    Info.stream().incrIndent();

    bool normal = dict.getOrDefault<bool>("normal","false");
    bool ok=true;
    

    if(normal)
    {
        scalar speed;
        ok &=dict.readEntry("velocity",speed);
        //Positive speed towards the disk
        velocity = - speed * rDiscrete->geometry().direction();
        indent(Info)<< "- Normal to rotor speed: "<<speed<<endl;
    }
    else
    {
        ok &=dict.readEntry("velocity",velocity);
    }

    indent(Info)<< "- Rotor input velocity: "<<velocity<<endl;

    forAll(this->sampledVel,i)
    {
        this->sampledVel[i]=velocity;
    }

    Info.stream().decrIndent();
    return ok;
    
}

const vectorField& fixedVelocity::sampleVelocity(const volVectorField& U)
{
    return this->sampledVel;
}

}
