#include "velocitySampler.H"
#include "runTimeSelectionTables.H"

namespace Foam
{
 
    defineTypeNameAndDebug(velocitySampler,0);
    defineRunTimeSelectionTable(velocitySampler, dictionary);

    void velocitySampler::writeSampled(const word& name)
    {
        volScalarField sampled
            (
                IOobject
                (
                    name + ":sampledCells",
                    rMesh->mesh().time().timeName(),
                    rMesh->mesh()
                ),
                rMesh->mesh(),
                dimensionedScalar(dimless, Zero)
        );
        sampled.write();
    }

    autoPtr<velocitySampler> velocitySampler::New(const dictionary &dict, const rotorDiscrete *rDiscrete_, const rotorFvMeshSel *rMesh_)
    {
        //Get model Type name (Ex: fixedVelocity) 
        //From type key from dictionary (propellerModel)
        const word modelType(dict.get<word>("type")); 
        Info<<endl;
        Info<< "Selecting " << typeName << " " << modelType << endl;

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

    velocitySampler::velocitySampler(const rotorDiscrete* rDiscrete_,const rotorFvMeshSel* rMesh_)
        : rDiscrete(rDiscrete_), rMesh(rMesh_), sampledVel(rDiscrete->grid.centers().size(),{0,0,0})
    {
        
    }




}

