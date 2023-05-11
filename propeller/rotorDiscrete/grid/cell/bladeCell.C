#include "bladeCell.H"


namespace Foam
{
bladeCell::bladeCell(scalar radius0, scalar radius1, scalar chord0, scalar chord1)
    : radius0_(radius0), radius1_(radius1), chord0_(chord0), chord1_(chord1)
{
    dr_ = this->radius1()-this->radius0();
    dc_ = this->chord1()-this->chord0();

    center_ = vector(this->radius0()+dr_/2,0,0);

    localPoints_.resize(4);

    localPoints_[0] = vector(this->radius0(),this->chord0(),0);
    localPoints_[1] = vector(this->radius1(),this->chord0(),0);    
    localPoints_[2] = vector(this->radius1(),this->chord1(),0);    
    localPoints_[3] = vector(this->radius0(),this->chord1(),0);    

    actualLocation_ = localPoints_;

}

void bladeCell::build()
{
    scalar totalw=0.0;
    Info<<cellis_<<endl;
    if(cellis_.size()==0)
    {
        FatalErrorInFunction
            <<"Some bladeCells doesnt contain any fvCell"
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

void bladeCell::setRotation(const tensor &rotation)
{
    forAll(localPoints_,i)
    {
        actualLocation_[i] = transform(rotation,localPoints_[i]);
    }
}


vector bladeCell::cartesianCenter() const
{
    vector c;
    forAll(actualLocation_,i)
    {
        c+=actualLocation_[i];
    }

    return c/actualLocation_.size();
}

void bladeCell::centerFromClosestCell(const vectorField &cellCenters, const coordSystem::cartesian &localCartesian)
{
    /*scalar minDistSq = VGREAT;
    vector gridCC = vector(radius0()+dr()/2,theta0()+dt()/2,0);
    gridCC = localCartesian.globalPosition(gridCC);
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
    center_ = localCartesian.localPosition(center_);
    center_.z()=0;*/
}

void bladeCell::setLocalTensor(const tensor& localTensor)
{
    localBlade_=localTensor;
}

vector bladeCell::scaleForce(const vector &globalForce)
{
   return globalForce * this->dr();
}
void bladeCell::applySource(vectorField &source, const scalarField& cellVol, vector& scaledForce) const
{
    forAll(cellis_,k)
    {
        source[cellis_[k]] =  -(weights_[k] * scaledForce/cellVol[cellis_[k]]);
    }
}

}
