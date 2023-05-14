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
                    rMesh_->mesh().time().timeName(),
                    rMesh_->mesh()
                ),
                rMesh_->mesh(),
                dimensionedScalar(dimless, Zero)
        );
        sampled.write();
    }

    autoPtr<velocitySampler> velocitySampler::New(const dictionary &dict, const rotorGrid *rGrid, const rotorFvMeshSel *rMesh)
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

        return autoPtr<Foam::velocitySampler>(ctorPtr(dict,rGrid,rMesh));

    }

    velocitySampler::velocitySampler(const rotorGrid* rGrid,const rotorFvMeshSel* rMesh)
        : rGrid_(rGrid), rMesh_(rMesh), sampledVel(rGrid_->nCells(),Zero)
    {
        
    }




}

