#include "diskSampler.H"
#include "runTimeSelectionTables.H"

namespace Foam
{


defineTemplateTypeNameWithName(diskSampler<scalar>,"scalarDiskSampler");
defineTemplatedRunTimeSelectionTable(diskSampler, dictionary, scalar); 


defineTemplateTypeNameWithName(diskSampler<vector>,"vectorDiskSampler");
defineTemplatedRunTimeSelectionTable(diskSampler, dictionary, vector); 

template<class fType>
void diskSampler<fType>::writeSampled(const word& name)
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
template<class fType>
autoPtr<diskSampler<fType>> diskSampler<fType>::New(const dictionary &dict, const rotorGrid *rGrid, const rotorFvMeshSel *rMesh)
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

    return autoPtr<Foam::diskSampler>(ctorPtr(dict,rGrid,rMesh));

}

template<class fType>
diskSampler<fType>::diskSampler(const rotorGrid* rGrid,const rotorFvMeshSel* rMesh)
    : rGrid_(rGrid), rMesh_(rMesh), sampledVel(rGrid_->nCells(),Zero)
{
    
}




}

