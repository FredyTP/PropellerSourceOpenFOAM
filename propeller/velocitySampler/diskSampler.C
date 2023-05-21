#include "diskSampler.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

//defineTemplateTypeNameAndDebugWithName(diskSampler<scalar>,"diskSampler",0);
//defineTemplatedRunTimeSelectionTable(diskSampler, dictionary, scalar); 

defineTemplateTypeNameAndDebugWithName(diskSampler<scalar>,"diskSampler",0);
defineTemplatedRunTimeSelectionTable(diskSampler, dictionary, scalar); 

defineTemplateTypeNameAndDebugWithName(diskSampler<vector>,"diskSampler",0);
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
        Info<<"No ptr"<<endl;
        FatalIOErrorInLookup
        (
            dict,
            typeName,
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    Info<<"Returning disk sampler"<<endl;
    return autoPtr<Foam::diskSampler<fType>>(ctorPtr(dict,rGrid,rMesh));

}

template<class fType>
diskSampler<fType>::diskSampler(const rotorGrid* rGrid,const rotorFvMeshSel* rMesh)
    : rGrid_(rGrid), rMesh_(rMesh), sampledField_(rGrid_->nCells(),Zero)
{
    
}


}