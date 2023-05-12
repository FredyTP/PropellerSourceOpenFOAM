#include "bladeCell.H"


namespace Foam
{
bladeCell::bladeCell(const rotorGeometry& rotorGeometry, scalar radius0, scalar radius1, scalar chord0, scalar chord1)
    : gridCell(rotorGeometry), radius0_(radius0), radius1_(radius1), chord0_(chord0), chord1_(chord1)
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

    factor_ = dr_;

}

bladeCell::bladeCell(const rotorGeometry &rotorGeometry, const List<vector> &localVertex)
    : gridCell(rotorGeometry)
{
    vector maxV = max(localVertex);
    vector minV = min(localVertex);
    dr_ = maxV.x()-minV.x();
    dc_=0;
    factor_ = dr_;
    localPoints_= localVertex;
    actualLocation_ = localPoints_;

}

void bladeCell::setRotation(const tensor &rotation)
{
    forAll(localPoints_,i)
    {
        actualLocation_[i] = transform(rotation,localPoints_[i]);
    }
    //Clear selection after rotation
    cellis_.clear();
    weights_.clear();
}


vector bladeCell::getCellCenter() const
{
    vector c = Zero;
    forAll(actualLocation_,i)
    {
        c+=actualLocation_[i];
    }

    return rotorGeometry_.cartesianToCylindrical().localPosition(c/actualLocation_.size());
}


}
