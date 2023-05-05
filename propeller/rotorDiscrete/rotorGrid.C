#include "rotorGrid.H"
#include "mathematicalConstants.H"
namespace Foam
{
rotorGrid::rotorGrid(label nRadius, label nTheta, scalar minRadius, scalar maxRadius)
    : minRadius_(minRadius),maxRadius_(maxRadius),radius_(nRadius+1),theta_(nTheta+1),ijkAddressing(nRadius,nTheta,0)
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

void rotorGrid::assignFvCells(const coordSystem::cylindrical &cylCS, const vectorField &cellCenter, const scalarField& weights, const labelList& cellis)
{
    forAll(cellis,i)
    {
        labell celli = cellis[i];
        vector cc = cellCenter[celli];
        vector polar = cylCS.localPosition(cc);
        label ir=0;
        label unused;
        label it=0;
        int result = regularInterpolation<scalar,scalar,1>::FindIndex(polar.x(),this->radius(),ir,unused);
        regularInterpolation<scalar,scalar,1>::FindIndex(polar.y(),this->theta(),it,unused);

        if(result == 1)
        {
            this->cell(ir,it).addCelli(rotorCells_[i].celli(),weights[rotorCells_[i].celli()]);
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
    }
}

void rotorGrid::setCenterFromClosestCell(coordSystem::cylindrical &cylCS, const vectorField &cellCenter)
{
    forAll(cells_,i)
    {
        cells_[i].centerFromClosestCell(cellCenter,cylCS);
    }
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