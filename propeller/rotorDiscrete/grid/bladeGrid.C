#include "bladeGrid.H"
#include "unitConversion.H"
#include "axisAngleRotation.H"
#include "bladeCell.H"
#include "delaunayTriangulation.H"
namespace Foam
{
bladeGrid::bladeGrid(const rotorGeometry& geometry,const rotorFvMeshSel& rotorFvMeshSel, scalar chord, label nBlades, label nRadius, label nChord)
: 
    rotorGrid(geometry,rotorFvMeshSel,nBlades), 
    ijkAddressing(nBlades,nRadius,nChord), 
    chord_(chord),
    nRadius_(nRadius),
    nChord_(nChord),
    theta_(nBlades)
{
    cells_.resize(this->size());
    buildBlades();
    updateTheta(0);
    rotateBlades();
    assignFvCells();
    build();
}
void bladeGrid::assignFvCells()
{
    const auto& cellCenter = meshSel_.mesh().C();
    const auto& weights = meshSel_.mesh().V();
    const auto& cellis = meshSel_.cells();

    const auto& carCS = rotorGeometry_.cartesianCS();
    forAll(cellis,i)
    {
        label celli = cellis[i];
        vector center = carCS.localPosition(cellCenter[celli]);
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
    const auto& cylCS = rotorGeometry_.cylindricalCS();
    centers_.resize(cells_.size());
    forAll(cells_,i)
    {
        cells_[i].checkCells();
        cells_[i].buildWeigths();
        vector tmp = cells_[i].getCellCenter();
        cells_[i].setCenter(tmp);
        centers_[i]=cells_[i].center();
    }
}
void bladeGrid::setRotation(scalar theta0)
{
    updateTheta(theta0);
    rotateBlades();
    assignFvCells();
    build();

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
    for(label ib = 0; ib < nBlades_ ;ib ++ )
    {
        for(label ir=0; ir<nRadius_; ir++)
        {
            for(label ic = 0; ic < nChord_; ic++)
            {
                bladeCell* newcell = new bladeCell(rotorGeometry_,radius[ir],radius[ir+1],chord[ic],chord[ic+1]);
                cells_.set(index(ib,ir,ic),newcell);
            }
        }
    }


}
void bladeGrid::updateTheta(scalar theta0)
{
    scalar dtheta = constant::mathematical::twoPi/nBlades_;
    forAll(theta_,i)
    {
        theta_[i]=theta0 + dtheta*i;

        if(theta_[i]>constant::mathematical::pi)
        {
            theta_[i]-=constant::mathematical::twoPi;
        }
    }
}

void bladeGrid::rotateBlades()
{
    for(label ib = 0; ib < nBlades_ ;ib ++ )
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
                    //vector center = bCell->cartesianCenter();
                    //center = rotorGeometry_.cartesianToCylindrical().localPosition(center);
                    //bCell->setCenter(center);
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
