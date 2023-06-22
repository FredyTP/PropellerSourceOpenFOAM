#include "azimutalMean.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolationCellPoint.H"

namespace Foam
{
 
defineTemplateTypeNameWithName(azimutalMean<scalar>,"azimutalMean");
addTemplatedToRunTimeSelectionTable(diskSampler,azimutalMean,scalar,dictionary);

defineTemplateTypeNameWithName(azimutalMean<vector>,"azimutalMean");
addTemplatedToRunTimeSelectionTable(diskSampler,azimutalMean,vector,dictionary);


template<class fType>
azimutalMean<fType>::azimutalMean(const dictionary& dict,const rotorGrid* rGrid,const rotorFvMeshSel* rMesh)
 : domainSampler<fType>(dict,rGrid,rMesh)
{
    polarGrid_ = dynamic_cast<const polarGrid*>(rGrid);
    if(polarGrid_==nullptr)
    {
        FatalErrorInFunction<< "azimutalMean velocitySampler is only available for rotorGrid: polarGrid"
        <<exit(FatalError);
    }
}

template<class fType>
const Field<fType>& azimutalMean<fType>::sampleField(const GeometricField<fType, fvPatchField, volMesh>& U) 
{
    domainSampler<fType>::sampleField(U);
    
    label nR = polarGrid_->nRadial();
    label nA = polarGrid_->nAzimutal();

    for(label iR = 0; iR < nR; iR++)
    {
        fType umean = Zero;
        for(label iA = 0; iA < nA; iA++)
        {
            label idx = polarGrid_->index(iR,iA,0);
            umean += this->sampledField_[idx];
        }
        umean/=nA;
        for(label iA = 0; iA < nA; iA++)
        {
            label idx = polarGrid_->index(iR,iA,0);
            this->sampledField_[idx] = umean;
        }
    }


    return this->sampledField_;

}


}
