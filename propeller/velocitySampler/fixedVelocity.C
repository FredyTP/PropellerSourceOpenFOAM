#include "fixedVelocity.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam
{
 
    defineTypeNameAndDebug(fixedVelocity,0);
    addToRunTimeSelectionTable(velocitySampler,fixedVelocity, dictionary);


fixedVelocity::fixedVelocity(const dictionary& dict,const rotorGrid* rGrid,const rotorFvMeshSel* rMesh)
    : velocitySampler(rGrid,rMesh)
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
        velocity = - speed * rGrid_->geometry().direction().get();
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
