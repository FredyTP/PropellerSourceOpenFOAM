#include "offsetSampler.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolationCellPoint.H"

namespace Foam
{
 
    defineTypeNameAndDebug(offsetSampler,0);
    addToRunTimeSelectionTable(velocitySampler,offsetSampler, dictionary);

offsetSampler::offsetSampler(const dictionary& dict,const rotorDiscrete* rDiscrete_,const rotorMesh* rMesh_)
 : velocitySampler(rDiscrete_,rMesh_)
{
    this->read(dict);
}

bool offsetSampler::read(const dictionary &dict)
{
    offset = dict.getOrDefault<scalar>("offset",0.0);
    if(std::abs(offset)<=SMALL)
    {
        offset=0.0;
    }
    atCellCenter = dict.getOrDefault<bool>("atCellCenter",true);
    
    Info<< "offset = "<<offset<<endl;
    Info<< "sample atCellCenter = "<<atCellCenter<<endl;

    this->build();

    return true;
}
 
const vectorField& offsetSampler::sampleVelocity(const volVectorField& U) 
{
    //If no offset and rotorDiscrete is equal to rotorMesh
    //Then the correspondence is cell to cell
    if(offset == 0.0 && rDiscrete->mode() == rotorDiscrete::dmMesh)
    {
        forAll(this->sampledVel,i)
        {
            this->sampledVel[i]=U.primitiveField()[rMesh->cells()[i]];
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
        //TODO:
        //This can be further improved finding cellweights when building
        //Update:
        //Using cellweights still extreme slow because the constructor
        //Interpolate on every cell face
        //This function should return a reference to a vectorField
        interpolationCellPoint<vector> interp(U);
        forAll(this->sampledVel,i)
        {
            this->sampledVel[i]=interp.interpolate(*(cellWeights[i].get()));
        }        
    }

    return this->sampledVel;
}
bool offsetSampler::build()
{
    //If offset is 0.0 and rotorDiscrete is equal to rotorMesh
    //There is no need to find cells or offset position, and the returned
    //velocity will be the velocity at cell center i of the rotor
    if(offset == 0.0 && rDiscrete->mode() == rotorDiscrete::dmMesh)
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
        //Get global coordinates
        point rPoint = rDiscrete->cylindrical().globalPosition(cylPoints[i]);

        //Add the offset normal to the geometry
        rPoint += rDiscrete->geometry().direction * offset;
        //Find the cell where the point is and set to the list
        cellToSample[i] = rMesh->mesh().findCell(rPoint); 

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
        sampled[cellToSample[i]]=1.0;
    }
    sampled.write();
}
}

