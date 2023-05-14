#include "domainSampler.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolationCellPoint.H"

namespace Foam
{
 
    defineTypeNameAndDebug(domainSampler,0);
    addToRunTimeSelectionTable(velocitySampler,domainSampler, dictionary);

domainSampler::domainSampler(const dictionary& dict,const rotorGrid* rGrid,const rotorFvMeshSel* rMesh)
 : velocitySampler(rGrid,rMesh)
{
    this->read(dict);
}

bool domainSampler::read(const dictionary &dict)
{
    offset = dict.getOrDefault<scalar>("offset",0.0);
    scale = dict.getOrDefault<scalar>("scale",1.0);
    if(std::abs(offset)<=VSMALL)
    {
        offset=0.0;
    }
    atCellCenter = dict.getOrDefault<bool>("atCellCenter",rGrid_->samplingMode() == rotorGrid::sampleMode::spClosestCell);
    Info.stream().incrIndent();
    Info<<indent<< "- Offset: "<<offset<<endl;
    Info<<indent<< "- Scale: "<<scale<<endl;
    Info<<indent<< "- Sample atCellCenter: "<<atCellCenter<<endl;

    //For parallel computation only 0 offset anc cell-center integration is available
    if(Pstream::parRun())
    {
        if(offset != 0.0)
        {
            Info<<indent<<"In parallel runs, only 0 offset and cell center integration is available"<<endl;
            FatalErrorInFunction<<exit(FatalError);
        }
    }
    
    this->build();

    Info.stream().decrIndent();
    return true;
}
 
const vectorField& domainSampler::sampleVelocity(const volVectorField& U) 
{
    //If no offset and rotorGrid is integrated in cell centers
    //Then the correspondence is cell to cell
    if(isDirectSample())
    {
        const PtrList<gridCell>& rotorCells = rGrid_->cells();
        forAll(rotorCells,i)
        {
            this->sampledVel[i] = U.primitiveField()[rotorCells[i].interpolatingCelli()];
        } 
    }
    else if(atCellCenter)
    {
        forAll(this->sampledVel,i)
        {
            this->sampledVel[i]=U.primitiveField()[cellToSample[i]];
        }     
    }
    else
    {
        interpolationCellPoint<vector> interp(U);
        forAll(this->sampledVel,i)
        {
            this->sampledVel[i]=interp.interpolate(*(cellWeights[i].get()));       
        }      
    }

    return this->sampledVel;
}
bool domainSampler::build()
{
    //If offset is 0.0 and rotorGrid is equal to rotorFvMeshSel
    //There is no need to find cells or offset position, and the returned
    //velocity will be the velocity at cell center i of the rotor
    Info<<"directSampler: "<<isDirectSample()<<endl;
    
    if(isDirectSample())
    {
        return true;
    }
    const auto& cells = rGrid_->cells();
    cellToSample.resize(cells.size());
    if(!atCellCenter)
    {
        cellWeights.resize(cells.size());
    }
    //Iterate over all discretization points
    forAll(cells, i)
    {
        //Get global coordinates
        vector localPoint = cells[i].center();
        //Scale radius
        localPoint.x() = localPoint.x()*scale;
        point rPoint = rGrid_->geometry().cylindricalCS().globalPosition(localPoint);
        
        //Add the offset normal to the geometry
        rPoint += rGrid_->geometry().direction().get() * offset;
        //Find the cell where the point is and set to the list
        cellToSample[i] = rMesh_->mesh().findCell(rPoint); 

        if(cellToSample[i]==-1)
        {
            FatalErrorInFunction
                << "sampledCell at position = "
                << localPoint 
                <<", is outside computational domain"
                <<exit(FatalError);
        }

        if(!atCellCenter)
        {
            cellWeights[i] = autoPtr<cellPointWeight>::New(rMesh_->mesh(),rPoint,cellToSample[i]);
        }
    }

    return true;
}

void domainSampler::writeSampled(const word& name)
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
    if(isDirectSample())
    {
        forAll(rGrid_->cells(),i)
        {
            const auto& cell = rGrid_->cells()[i];
            sampled[cell.interpolatingCelli()] +=1.0;
        }
    }
    else
    {
        forAll(cellToSample,i)
        {
            sampled[cellToSample[i]]+=1.0;   
        }
    }

    sampled.write();
}
bool domainSampler::isDirectSample()
{
    return (offset == 0.0 
    && scale == 1.0
    && atCellCenter 
    && rGrid_->samplingMode() == rotorGrid::sampleMode::spClosestCell);
}

}
