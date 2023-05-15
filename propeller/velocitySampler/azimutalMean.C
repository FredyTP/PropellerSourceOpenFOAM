#include "azimutalMean.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolationCellPoint.H"

namespace Foam
{
 
    defineTypeNameAndDebug(azimutalMean,0);
    addToRunTimeSelectionTable(velocitySampler,azimutalMean, dictionary);

azimutalMean::azimutalMean(const dictionary& dict,const rotorGrid* rGrid,const rotorFvMeshSel* rMesh)
 : domainSampler(dict,rGrid,rMesh)
{
    polarGrid_ = dynamic_cast<const polarGrid*>(rGrid);
    if(polarGrid_==nullptr)
    {
        FatalErrorInFunction<< "azimutalMean velocitySampler is only available for rotorGrid: polarGrid"
        <<exit(FatalError);
    }
}
 
const vectorField& azimutalMean::sampleVelocity(const volVectorField& U)
{
    domainSampler::sampleVelocity(U);
    
    label nR = polarGrid_->nRadial();
    label nA = polarGrid_->nAzimutal();

    for(label iR = 0; iR < nR; iR++)
    {
        vector umean = Zero;
        for(label iA = 0; iA < nA; iA++)
        {
            label idx = polarGrid_->index(iR,iA,0);
            umean += this->sampledVel[idx];
        }
        umean/=nA;
        for(label iA = 0; iA < nA; iA++)
        {
            label idx = polarGrid_->index(iR,iA,0);
            this->sampledVel[idx] = umean;
        }
    }

    if(rMesh_->mesh().time().writeTime())
    {
        volVectorField vel
        (
            IOobject
            (
                "sampledVel",
                rMesh_->mesh().time().timeName(),
                rMesh_->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rMesh_->mesh(),
            dimensionedVector(dimVelocity, Zero)
        );
        forAll(polarGrid_->cells(),i)
        {
            polarGrid_->cells()[i].applyField<vector>(vel,sampledVel[i]);
        }
        vel.write();
    }
    
    return sampledVel;

}


}
