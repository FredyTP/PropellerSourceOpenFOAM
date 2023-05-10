#include "bladeGrid.H"
#include "unitConversion.H"
#include "axisAngleRotation.H"
#include "bladeCell.H"
#include "delaunayTriangulation.H"
namespace Foam
{
bladeGrid::bladeGrid(const coordSystem::cartesian &carCS, const coordSystem::cylindrical &cylCS, scalar innerRadius, scalar radius,scalar chord, label nBlade, label nRadius, label nChord)
: 
    rotorGrid(cylCS,innerRadius,radius), 
    ijkAddressing(nBlade,nRadius,nChord), 
    carCS_(carCS),
    chord_(chord),
    nBlade_(nBlade), 
    nRadius_(nRadius),
    nChord_(nChord),
    theta_(nBlade)
{
    cells_.resize(this->size());
    buildBlades();
    updateTheta(0);
    rotateBlades();
}
void bladeGrid::assignFvCells(const vectorField &cellCenter, const scalarField &weights, const labelList &cellis)
{
    Info<<"assigning"<<endl;
    forAll(cells_,i)
    {
        if(cells_.get(i)==nullptr)
        {
            Info<<"cell: "<<i<<" is NULL"<<endl;
        }
    }
    forAll(cellis,i)
    {
        label celli = cellis[i];
        vector center = carCS_.localPosition(cellCenter[celli]);
        List<label> fakeCell({0,1,2,3});
        forAll(cells_,j)
        {
            bladeCell* bCell = dynamic_cast<bladeCell*>(cells_.get(j));
            if(bCell == nullptr)
            {
                continue;
            }       
            if(delaunayTriangulation::isInsideCell(bCell->actualPoints(),fakeCell,center))
            {
                cells_.get(j)->addCelli(celli,weights[celli]);            
                break;
            }
        }
    }
}
void bladeGrid::build()
{

}
void bladeGrid::buildBlades()
{
    List<scalar> radius(nRadius_+1);
    List<scalar> chord(nChord_+1);

    scalar dr = (maxRadius_-minRadius_)/nRadius_;
    scalar dc = (chord_)/nChord_;

    //Create radius points
    for(label i = 0; i <nRadius_+1;i++)
    {
        radius[i]=minRadius_+dr*i;
    }

    for(label i = 0; i <nChord_+1;i++)
    {
        chord[i]=-chord_/2 + dc*i;
    }

    cells_.resize(this->size());
    for(label ib = 0; ib < nBlade_ ;ib ++ )
    {
        for(label ir=0; ir<nRadius_; ir++)
        {
            for(label ic = 0; ic < nChord_; ic++)
            {
                bladeCell* newcell = new bladeCell(radius[ir],radius[ir+1],chord[ic],chord[ic+1]);
                cells_.set(index(ib,ir,ic),newcell);
                if(cells_.get(index(ib,ir,ic))==nullptr)
                {
                    Info<<"Cell: " <<index(ib,ir,ic)<<" is null"<<endl;
                }
                else
                {
                    Info<<"Cell: " <<index(ib,ir,ic)<<" is NOT null"<<endl;
                }
            }
        }
    }

}
void bladeGrid::updateTheta(scalar theta0)
{
    scalar dtheta = constant::mathematical::twoPi/nBlade_;
    forAll(theta_,i)
    {
        theta_[i]=theta0 + dtheta*i;

        if(theta_[i]>constant::mathematical::pi)
        {
            theta_[i]=theta_[i]-constant::mathematical::twoPi;
        }
    }
    Info<<"Thetas: "<<endl;
    Info<<theta_*radToDeg()<<endl;
}

void bladeGrid::rotateBlades()
{
    for(label ib = 0; ib < nBlade_ ;ib ++ )
    {
        tensor rotation = coordinateRotations::axisAngle::rotation(vector::components::Z,theta_[ib],false);
        for(label ir=0; ir<nRadius_; ir++)
        {
            for(label ic = 0; ic < nChord_; ic++)
            {
                bladeCell* bCell = dynamic_cast<bladeCell*>(cells_.get(index(ib,ir,ic)));
                if(bCell != nullptr)
                {
                    bCell->setRotation(rotation);
                }
            }
        }
    }
}
label bladeGrid::ijkIndex(label iBlade, label iRadius, label iChord)
{
    return ijkAddressing::index(iBlade,iRadius,iChord);
}

}
