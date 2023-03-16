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
    offset = dict.getOrDefault<bool>("offset",0);
    if(abs(offset)<=SMALL)
    {
        offset=0.0;
    }
    atCellCenter = dict.getOrDefault<bool>("atCellCenter",true);
    
    Info<< "offset = "<<offset<<endl;
    Info<< "sample atCellCenter = "<<atCellCenter<<endl;

    this->build();
    
    return true;
}
 
vector offsetSampler::sampleVelocityAt(const volVectorField &U, label i) const
{
    //If no offset and rotorDiscrete is equal to rotorMesh
    //Then the correspondence is cell to cell
    if(offset == 0.0 && rDiscrete->mode() == rotorDiscrete::dmMesh)
    {
        return U.primitiveField()[rMesh->cells()[i]];
    }
    else if(atCellCenter)
    {
        return U.primitiveField()[cellToSample[i]];
    }
    else
    {
        //TODO:
        //This can be further improved finding cellweights when building
        interpolationCellPoint<vector> interp(U);
        return interp.interpolate(posToSample[i],cellToSample[i]);
    }
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
    posToSample.resize(cylPoints.size());

    //Iterate over all discretization points
    forAll(cylPoints, i)
    {
        //Get global coordinates
        point rPoint = rDiscrete->cylindrical().globalPosition(cylPoints[i]);

        //Add the offset normal to the geometry
        rPoint += rDiscrete->geometry().direction * offset;
        posToSample[i] = rPoint;
        //Find the cell where the point is and set to the list
        cellToSample[i] = rMesh->mesh().findCell(rPoint); 
    }

    return true;
}
}

