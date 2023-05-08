#include "rotorGrid.H"
#include "mathematicalConstants.H"
#include "regularInterpolation.H"

namespace Foam
{
rotorGrid::rotorGrid(label nRadius, label nTheta, scalar minRadius, scalar maxRadius,const coordSystem::cylindrical &cylCS)
    : cylCS_(cylCS), minRadius_(minRadius),maxRadius_(maxRadius),radius_(nRadius+1),theta_(nTheta+1),ijkAddressing(nRadius,nTheta,0)
{
    scalar dr = (maxRadius_-minRadius_)/nRadius;
    scalar dt = (constant::mathematical::twoPi)/nTheta;


    for(label i = 0; i <nRadius+1;i++)
    {
        radius_[i]=minRadius_+dr*i;
    }


    for(label i = 0; i <nTheta+1;i++)
    {
        theta_[i]=-constant::mathematical::pi+dt*i;
    }
    //Closed loop

    buildGrid();
}

void rotorGrid::assignFvCells(const vectorField &cellCenter, const scalarField& weights, const labelList& cellis)
{
    forAll(cellis,i)
    {
        label celli = cellis[i];
        vector cc = cellCenter[celli];
        vector polar = cylCS_.localPosition(cc);
        label ir=0;
        label unused;
        label it=0;
        int result = regularInterpolation<scalar,scalar,1>::FindIndex(polar.x(),this->radius(),ir,unused);
        regularInterpolation<scalar,scalar,1>::FindIndex(polar.y(),this->theta(),it,unused);

        if(result == 1)
        {
            this->cell(ir,it).addCelli(celli,weights[celli]);
        }
    }
}

void rotorGrid::build()
{
    centers_.resize(cells_.size());
    forAll(cells_,i)
    {
        cells_[i].build();
        centers_[i]=cells_[i].center();
        tensor bladetensor = rotorGrid::bladeLocalFromPoint(cylCS_,cells_[i].center());
        cells_[i].setLocalTensor(bladetensor);
    }
}

void rotorGrid::setCenterFromClosestCell(const vectorField &cellCenter)
{
    forAll(cells_,i)
    {
        cells_[i].centerFromClosestCell(cellCenter,cylCS_);
    }
}

tensor rotorGrid::bladeLocalFromPoint(const coordSystem::cylindrical &cylCS, const point &localPoint) 
{
    // z- up, y -outwards from center, x perpendicular y,z (leading edge to trailing edge)
    point global, origin;

    global = cylCS.globalPosition(localPoint);
    origin = cylCS.origin();

    tensor rotTensor(cylCS.R());

    // z-axis
    const vector ax3 = rotTensor.col<2>(); // == e3 (already normalized)

    // y-axis (radial direction)
    vector ax2(global - origin);

    ax2.removeCollinear(ax3);

    const scalar magAxis2(mag(ax2));

    // Trap zero size and colinearity
    if (magAxis2 < SMALL)
    {
        return rotTensor;
    }

    ax2 /= magAxis2; // normalise

    // Replace with updated local axes

    rotTensor.col<0>(ax2 ^ ax3);
    rotTensor.col<1>(ax2);

    return rotTensor;
}
void rotorGrid::buildGrid()
{
    label radialCells = radius_.size()-1;
    label thetaCells = theta_.size()-1;

    cells_.resize(radialCells*thetaCells);

    for(label i = 0; i< radialCells;i++)
    {
        for(label j = 0; j<thetaCells;j++)
        {
            cells_.emplace(index(i,j,0),radius_[i],radius_[i+1],theta_[i],theta_[i+1]);
        }
    }
}

}