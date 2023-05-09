#include "gridCell.H"


namespace Foam
{
gridCell::gridCell(scalar radius0, scalar radius1, scalar theta0, scalar theta1)
    : radius0_(radius0), radius1_(radius1), theta0_(theta0), theta1_(theta1)
{
    dr_ = this->radius1()-this->radius0();
    dt_ = this->theta1()-this->theta0();

    center_ = vector(this->radius0()+dr_/2,this->theta0()+dt_/2,0);
}

void gridCell::addCelli(label celli,scalar weight)
{
    cellis_.append(celli);
    weights_.append(weight);
}
void gridCell::build()
{
    scalar totalw=0.0;
    if(cellis_.size()==0)
    {
        FatalErrorInFunction
            <<"Some gridCells doesnt contain any fvCell"
            <<exit(FatalError);
    }
    
    forAll(weights_,i)
    {
        totalw+=weights_[i];
    }

    forAll(weights_,i)
    {
        weights_[i]/=totalw;
    }
}

void gridCell::setCenter(vector &center)
{
    center_=center;
    center_.z()=0;
}

void gridCell::centerFromClosestCell(const vectorField &cellCenters, const coordSystem::cylindrical &localCyl)
{
    scalar minDistSq = VGREAT;
    vector gridCC = vector(radius0()+dr()/2,theta0()+dt()/2,0);
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

    center_ = cellCenters[interpolatingCell_];
    center_ = localCyl.localPosition(center_);
    center_.z()=0;
}

void gridCell::setLocalTensor(const tensor& localTensor)
{
    localBlade_=localTensor;
}

vector gridCell::scaleForce(const vector &globalForce)
{
   return globalForce * this->dr() * this->dt() / constant::mathematical::twoPi;
}
void gridCell::applySource(vectorField &source, const scalarField& cellVol, vector& scaledForce) const
{
    forAll(cellis_,k)
    {
        source[cellis_[k]] =  -(weights_[k] * scaledForce/cellVol[cellis_[k]]);
    }
}

}
