#include "fixedVelocity.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam
{
 
    defineTypeNameAndDebug(fixedVelocity,0);
    addToRunTimeSelectionTable(velocitySampler,fixedVelocity, dictionary);


fixedVelocity::fixedVelocity(const dictionary& dict,const rotorDiscrete* rDiscrete_,const rotorMesh* rMesh_)
    : velocitySampler(rDiscrete_,rMesh_)
{
    this->read(dict);
}

bool fixedVelocity::read(const dictionary &dict)
{
    bool normal = dict.getOrDefault<bool>("normal","false");
    bool ok=true;
    

    if(normal)
    {
        scalar speed;
        ok &=dict.readEntry("velocity",speed);
        //Positive speed towards the disk
        velocity = - speed * rDiscrete->geometry().direction;
        Info<< "normal to rotor speed = "<<speed<<endl;
    }
    else
    {
        ok &=dict.readEntry("velocity",velocity);
    }

    Info<< "rotor input velocity = "<<velocity<<endl;
    return ok;
    
}

vector fixedVelocity::sampleVelocityAt(const volVectorField &U, label i) const
{
    return velocity;
}
}
