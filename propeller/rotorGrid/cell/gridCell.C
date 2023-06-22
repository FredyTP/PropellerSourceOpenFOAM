#include "gridCell.H"
#include "rotorGrid.H"
#include "Pstream.H"
namespace Foam
{


void gridCell::addCelli(label celli,scalar weight)
{
    cellis_.append(celli);
    weights_.append(weight);
}


void gridCell::buildWeigths(bool parCheck)
{
    scalar totalw=0.0;

    forAll(weights_,i)
    {
        totalw+=weights_[i];
    }

    reduce(totalw,sumOp<scalar>());

    forAll(weights_,i)
    {
        weights_[i]/=totalw;
    }
}
void gridCell::checkCells(bool parCheck, label core)
{
    label ncell = cellis_.size();
    
    if(parCheck)
    {
        reduce(ncell,sumOp<label>());
        if(ncell==0)
        {
            FatalErrorInFunction
                <<"Some gridCells doesnt contain any fvCell"
                <<exit(FatalError);
        }
    }
    else if(core == Pstream::myProcNo())
    {
        if(ncell==0)
        {
            FatalErrorInFunction
                <<"Some gridCells doesnt contain any fvCell"
                <<exit(FatalError);
        }
    }

}


void gridCell::setCenter(const vector &center)
{
    center_=center;
    center_.z()=0;
}


void gridCell::centerFromClosestCell(const vectorField &cellCenters)
{
    const auto& localCyl = rotorGeometry_.cylindricalCS();
    scalar minDistSq = VGREAT;
    vector gridCC = getCellCenter();
    gridCC = localCyl.globalPosition(gridCC);
    interpolatingCell_ = -1;
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

    vector newCenter;
    Tuple2<scalar,vector> distpos;
    if(interpolatingCell_!=-1)
    {
        newCenter = cellCenters[interpolatingCell_];
        newCenter = localCyl.localPosition(newCenter);
        newCenter.z()=0;
    }

    distpos.first()=minDistSq;
    distpos.second()=newCenter;

    reduce(distpos,minFirstOp<scalar>());    
    this->setCenter(distpos.second());

    //Only interpolate cell from the closest one and disable other cores
    if(newCenter != distpos.second())
    {
        interpolatingCell_ = -1;
    }
    
}


vector gridCell::scaleForce(const vector &globalForce) const
{
   return globalForce * factor_;
}


}
