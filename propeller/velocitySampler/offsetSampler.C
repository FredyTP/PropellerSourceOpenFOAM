#include "offsetSampler.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolationCellPoint.H"

namespace Foam
{
 
    defineTypeNameAndDebug(offsetSampler,0);
    addToRunTimeSelectionTable(velocitySampler,offsetSampler, dictionary);

offsetSampler::offsetSampler(const dictionary& dict,const rotorDiscrete* rDiscrete_,const rotorFvMeshSel* rMesh_)
 : velocitySampler(rDiscrete_,rMesh_)
{
    this->read(dict);
}

bool offsetSampler::read(const dictionary &dict)
{
    offset = dict.getOrDefault<scalar>("offset",0.0);
    if(std::abs(offset)<=VSMALL)
    {
        offset=0.0;
    }
    atCellCenter = dict.getOrDefault<bool>("atCellCenter",true);
    Info.stream().incrIndent();
    Info<<indent<< "- Offset: "<<offset<<endl;
    Info<<indent<< "- Sample atCellCenter: "<<atCellCenter<<endl;

    //For parallel computation only 0 offset anc cell-center integration is available
    if(Pstream::parRun())
    {
        if(offset != 0.0 || rDiscrete->integrationMode() != rotorCell::integrationMode::imCenter)
        {
            Info<<indent<<"In parallel runs, only 0 offset and cell center integration is available"<<endl;
            FatalErrorInFunction<<exit(FatalError);
        }
    }
    
    this->build();

    Info.stream().decrIndent();
    return true;
}
 
const vectorField& offsetSampler::sampleVelocity(const volVectorField& U) 
{
    //If no offset and rotorDiscrete is integrated in cell centers
    //Then the correspondence is cell to cell
    if(offset == 0.0 && rDiscrete->integrationMode() == rotorCell::integrationMode::imCenter)
    {
        const PtrList<rotorCell>& rotorCells = rDiscrete->rotorCells();
        forAll(rotorCells,i)
        {
            this->sampledVel[rotorCells[i].center()] = U.primitiveField()[rotorCells[i].celli()];
        } 
    }
    else if(atCellCenter)
    {
        forAll(this->sampledVel,i)
        {
            if(rDiscrete->integrationPoints()[i])
            {
                this->sampledVel[i]=U.primitiveField()[cellToSample[i]];
            }
        }     
    }
    else
    {
        interpolationCellPoint<vector> interp(U);
        forAll(this->sampledVel,i)
        {
            if(rDiscrete->integrationPoints()[i])
            {
                this->sampledVel[i]=interp.interpolate(*(cellWeights[i].get()));
            }          
        }        
    }

    return this->sampledVel;
}
bool offsetSampler::build()
{
    //If offset is 0.0 and rotorDiscrete is equal to rotorFvMeshSel
    //There is no need to find cells or offset position, and the returned
    //velocity will be the velocity at cell center i of the rotor
    if(offset == 0.0 && rDiscrete->integrationMode() == rotorCell::integrationMode::imCenter)
    {
        return true;
    }
    const List<point>& cylPoints = rDiscrete->cylPoints();
    cellToSample.resize(cylPoints.size());
    if(!atCellCenter)
    {
        cellWeights.resize(cylPoints.size());
    }


    //Iterate over all discretization points
    forAll(cylPoints, i)
    {
        //Just get on the integration points
        if(!rDiscrete->integrationPoints()[i])
        {
            continue;
        }
        //Get global coordinates
        point rPoint = rDiscrete->cylindrical().globalPosition(cylPoints[i]);

        //Add the offset normal to the geometry
        rPoint += rDiscrete->geometry().direction() * offset;
        //Find the cell where the point is and set to the list
        cellToSample[i] = rMesh->mesh().findCell(rPoint); 

        if(cellToSample[i]==-1)
        {
            FatalErrorInFunction
                << "sampledCell at position = "
                << rPoint 
                <<", is outside computational domain"
                <<exit(FatalError);
                
        }

        if(!atCellCenter)
        {
            cellWeights[i] = autoPtr<cellPointWeight>::New(rMesh->mesh(),rPoint,cellToSample[i]);
        }
    }

    return true;
}

void offsetSampler::writeSampled(const word& name)
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

    forAll(cellToSample,i)
    {
        if(!rDiscrete->integrationPoints()[i])
        {
            continue;
        }

        sampled[cellToSample[i]]=1.0;   

        
    }
    sampled.write();
}
}

