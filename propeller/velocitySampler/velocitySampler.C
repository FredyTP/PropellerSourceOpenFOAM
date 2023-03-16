#include "velocitySampler.H"
#include "runTimeSelectionTables.H"

namespace Foam
{
 
    defineTypeNameAndDebug(velocitySampler,0);
    defineRunTimeSelectionTable(velocitySampler, dictionary);

    autoPtr<velocitySampler> velocitySampler::New(const dictionary &dict,const rotorDiscrete* rDiscrete_,const rotorMesh* rMesh_)
    {
        //Get model Type name (Ex: fixedVelocity) 
        //From type key from dictionary (propellerModel)
        const word modelType(dict.get<word>("type")); 

        Info<< "    Selecting " << typeName << " " << modelType << endl;

        //Find class contructor in tables
        auto* ctorPtr = dictionaryConstructorTable(modelType);

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                dict,
                typeName,
                modelType,
                *dictionaryConstructorTablePtr_
            ) << exit(FatalIOError);
        }

        return autoPtr<Foam::velocitySampler>(ctorPtr(dict,rDiscrete_,rMesh_));

    }

    velocitySampler::velocitySampler(const rotorDiscrete* rDiscrete_,const rotorMesh* rMesh_)
        : rDiscrete(rDiscrete_),rMesh(rMesh_)
    {
    }




}

