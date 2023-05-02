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

void rotorGrid::build()
{
    centers_.resize(cells_.size());
    forAll(cells_,i)
    {
        cells_[i].build();
        centers_[i]=cells_[i].center();
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
            cells_.emplace(index(i,j,0),i,j,radius_,theta_);
        }
    }
}

}