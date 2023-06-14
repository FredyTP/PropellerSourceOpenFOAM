#include "domainSampler.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolationCellPoint.H"

namespace Foam
{
 


defineTemplateTypeNameWithName(domainSampler<scalar>,"domainSampler");
addTemplatedToRunTimeSelectionTable(diskSampler,domainSampler,scalar,dictionary);

defineTemplateTypeNameWithName(domainSampler<vector>,"domainSampler");
addTemplatedToRunTimeSelectionTable(diskSampler,domainSampler,vector,dictionary);



template<class fType>
domainSampler<fType>::domainSampler(const dictionary& dict,const rotorGrid* rGrid,const rotorFvMeshSel* rMesh)
 : diskSampler<fType>(rGrid,rMesh)
{
    this->read(dict);
}
template<class fType>
bool domainSampler<fType>::read(const dictionary &dict)
{
    offset = dict.getOrDefault<scalar>("offset",0.0);
    scale = dict.getOrDefault<scalar>("scale",1.0);
    if(std::abs(offset)<=VSMALL)
    {
        offset=0.0;
    }
    Info.stream().incrIndent();
    Info<<indent<< "- Offset: "<<offset<<endl;
    Info<<indent<< "- Scale: "<<scale<<endl;
    Info<<indent<< "- Sample Mode: "<<rotorGrid::sampleModeNames_.get(this->rGrid_->samplingMode()) <<endl;
    Info<<indent<< "- Direct Sample: "<< isDirectSample() <<endl;

    if(!isDirectSample() && this->rGrid_->samplingMode() == rotorGrid::sampleMode::spCellMean)
    {
        FatalErrorInFunction<<"Cell Mean sampling mode is only available for 0.0 offset and 1.0 scale"
        <<exit(FatalError);
    }
    //For parallel computation only 0 offset and cell-center integration is available
    /*if(Pstream::parRun() && this->rMesh_->isMultiCore())
    {
        if(!isDirectSample())
        {
            Info<<indent<<"In parallel runs, only 0 offset is available"<<endl;
            //FatalErrorInFunction<<exit(FatalError);
        }
    }*/
    
    this->build();

    Info.stream().decrIndent();
    return true;
}
 
template<class fType>
const Field<fType>& domainSampler<fType>::sampleField(const GeometricField<fType, fvPatchField, volMesh>& U) 
{
    //If no offset and rotorGrid is integrated in cell centers
    //Then the correspondence is cell to cell
    if(isDirectSample())
    {
        const PtrList<gridCell>& rotorCells = this->rGrid_->cells();

        if(this->rGrid_->samplingMode() == rotorGrid::sampleMode::spClosestCell)
        {
            forAll(rotorCells,i)
            {
                label celli = rotorCells[i].interpolatingCelli();
                if(celli!=-1)
                {
                    this->sampledField_[i] = U.primitiveField()[celli];
                }
                else
                {
                    this->sampledField_[i] = Zero;
                }
            }       
        }
        else //Average
        {
            forAll(rotorCells,i)
            {
                this->sampledField_[i] = rotorCells[i].applyWeights<fType>(U.primitiveField());
            } 
            
        }

    }
    else if(this->rGrid_->samplingMode() == rotorGrid::sampleMode::spClosestCell)
    {
        forAll(this->sampledField_,i)
        {
            label celli = cellToSample[i];
            if(celli != -1)
            {
                this->sampledField_[i]=U.primitiveField()[celli];
            }
            else
            {
                this->sampledField_[i]=Zero;
            }
        }     

    }
    else
    {
        interpolationCellPoint<fType> interp(U);
        forAll(this->sampledField_,i)
        {
            if(cellWeights[i].good())
            {
                this->sampledField_[i]=interp.interpolate(*(cellWeights[i].get()));       
            }
            else
            {
                this->sampledField_[i] = Zero;
            }
        }      
    }
    reduce(this->sampledField_,sumOp<Field<fType>>());

    return this->sampledField_;
}
template<class fType>
void domainSampler<fType>::build()
{
    //If offset is 0.0 and rotorGrid is equal to rotorFvMeshSel
    //There is no need to find cells or offset position, and the returned
    //velocity will be the velocity at cell center i of the rotor
    if(isDirectSample())
    {
        return;
    }
    const auto& cells = this->rGrid_->cells();
    cellToSample.resize(cells.size());
    if(this->rGrid_->samplingMode() == rotorGrid::sampleMode::spCenter)
    {
        cellWeights.resize(cells.size());
    }

    //Check if cell found in any of the cores
    labelList foundCell(cells.size());
    //Iterate over all discretization points
    forAll(cells, i)
    {
        //Get global coordinates
        vector localPoint = cells[i].center();

        //Scale radius
        localPoint.x() = localPoint.x()*scale;
        point rPoint = this->rGrid_->geometry().cylindricalCS().globalPosition(localPoint);
        
        //Add the offset normal to the geometry
        rPoint += this->rGrid_->geometry().direction().get() * offset;
        
        //Find the cell where the point is and set to the list
        cellToSample[i] = this->rMesh_->mesh().findCell(rPoint); 

        if(cellToSample[i]==-1)
        {
            foundCell[i] = 0;
        }
        else
        {
            foundCell[i] = 1;
        }

        if(this->rGrid_->samplingMode() == rotorGrid::sampleMode::spCenter)
        {
            if(cellToSample[i] != -1)
            {
                cellWeights[i] = autoPtr<cellPointWeight>::New(this->rMesh_->mesh(),rPoint,cellToSample[i]);
            }
        }
    }
    
    //Check cellToSample
    reduce(foundCell,sumOp<labelList>());
    label idx = foundCell.find(0);
    if(idx != -1)
    {
        FatalErrorInFunction
        << "sampledCell at position = "
        << cells[idx].center()
        <<", is outside computational domain"
        <<exit(FatalError);
    }
    
}
template<class fType>
void domainSampler<fType>::writeSampled(const word& name)
{
    volScalarField sampled
        (
            IOobject
            (
                name + ":sampledCells",
                this->rMesh_->mesh().time().timeName(),
                this->rMesh_->mesh()
            ),
            this->rMesh_->mesh(),
            dimensionedScalar(dimless, Zero)
    );
    if(isDirectSample())
    {
        if(this->rGrid_->samplingMode() == rotorGrid::sampleMode::spClosestCell)
        {
            forAll(this->rGrid_->cells(),i)
            {
                const auto& cell = this->rGrid_->cells()[i];
                label celli = cell.interpolatingCelli();
                if(celli!=-1)
                {
                    sampled[celli] +=1.0;
                }
            }
        }
        else
        {
            forAll(this->rGrid_->cells(),i)
            {
                const auto& cell = this->rGrid_->cells()[i];
                forAll(cell.cellis(),j)
                {
                    sampled[cell.cellis()[j]] +=1.0;
                }
                
            }
        }

    }
    else
    {
        forAll(cellToSample,i)
        {
            label celli = cellToSample[i];
            if(celli != -1)
            {
                sampled[celli]+=1.0;   
            }
        }
    }

    sampled.write();
}
template<class fType>
bool domainSampler<fType>::isDirectSample()
{
    return (offset == 0.0 
    && scale == 1.0
    && this->rGrid_->samplingMode() != rotorGrid::sampleMode::spCenter);
}

}
