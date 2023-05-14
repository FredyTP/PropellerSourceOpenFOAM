#include "gridCell.H"
#include "rotorGrid.H"

namespace Foam
{


void gridCell::addCelli(label celli,scalar weight)
{
    cellis_.append(celli);
    weights_.append(weight);
}


void gridCell::buildWeigths()
{
    scalar totalw=0.0;

    forAll(weights_,i)
    {
        totalw+=weights_[i];
    }

    forAll(weights_,i)
    {
        weights_[i]/=totalw;
    }
}
void gridCell::checkCells()
{
    if(cellis_.size()==0)
    {
        FatalErrorInFunction
            <<"Some gridCells doesnt contain any fvCell"
            <<exit(FatalError);
    }
}


void gridCell::setCenter(const vector &center)
{
    center_=center;
    center_.z()=0;
    updateBladeTensor();
}

void gridCell::updateBladeTensor()
{
    localBlade_ = rotorGrid::bladeLocalFromPoint(rotorGeometry_.cylindricalCS(),center_);
}

void gridCell::centerFromClosestCell(const vectorField &cellCenters)
{
    const auto& localCyl = rotorGeometry_.cylindricalCS();
    scalar minDistSq = VGREAT;
    vector gridCC = getCellCenter();
    gridCC = localCyl.globalPosition(gridCC);
    forAll(cellis_,i)
    {
        label celli = cellis_[i];
        scalar distSq = magSqr(cellCenters[celli]-gridCC);
        
        if(distSq<minDistSq)
        {
            minDistSq=distSq;
            interpolatingCell_ = celli;
        }
    }
    vector newCenter = cellCenters[interpolatingCell_];
    newCenter = localCyl.localPosition(newCenter);
    newCenter.z()=0;
    this->setCenter(newCenter);
}


vector gridCell::scaleForce(const vector &globalForce)
{
   return globalForce * factor_;
}
void gridCell::applySource(vectorField &source, const scalarField& cellVol, vector& scaledForce) const
{
    forAll(cellis_,k)
    {
        source[cellis_[k]] =  -(weights_[k] * scaledForce/cellVol[cellis_[k]]);
    }
}

}
