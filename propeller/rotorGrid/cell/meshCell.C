#include "meshCell.H"
#include "geometry.H"

namespace Foam
{

meshCell::meshCell(const rotorGeometry& rotorGeometry, label celli, label nBlades, const List<vector>& localVertex, const vector* center)
: gridCell(rotorGeometry), localPoints_(localVertex)
{
    cellis_.append(celli);
    weights_.append(1.0);
    interpolatingCell_ = celli;
    if(center!=nullptr)
    {
        this->setCenter(rotorGeometry_.cartesianToCylindrical().localPosition(*center));
    }
    else
    {
        this->setCenter(getCellCenter());
    }

    util::geometry::sortCounterClockwise(localPoints_);
    area_ = util::geometry::poligonArea(localPoints_);
    factor_ = nBlades*area_/(constant::mathematical::twoPi * this->center().x() + SMALL);
}

meshCell::meshCell(const rotorGeometry &rotorGeometry, label celli, label nBlades, const List<point> &points, const List<label> &vertexIndex, const vector *center)
: gridCell(rotorGeometry)
{
    cellis_.append(celli);
    interpolatingCell_ = celli;
    weights_.append(1.0);

    localPoints_.resize(vertexIndex.size());
    forAll(localPoints_,i)
    {
        localPoints_[i] = points[vertexIndex[i]];
    }

    if(center!=nullptr)
    {
        this->setCenter(rotorGeometry_.cartesianToCylindrical().localPosition(*center));
    }
    else
    {
        this->setCenter(getCellCenter());
    }

    util::geometry::sortCounterClockwise(localPoints_);
   
    area_ = util::geometry::poligonArea(localPoints_);
    factor_ = nBlades*area_/(constant::mathematical::twoPi * this->center().x() + SMALL);
}
meshCell::meshCell(const rotorGeometry &rotorGeometry, label celli, label nBlades, scalar area, const vector &center)
: gridCell(rotorGeometry)
{
    cellis_.append(celli);
    interpolatingCell_ = celli;
    weights_.append(1.0);
    this->setCenter(rotorGeometry_.cartesianToCylindrical().localPosition(center));
    factor_ = nBlades*area_/(constant::mathematical::twoPi * this->center().x() + SMALL);
    area_ = area;
}

vector meshCell::getCellCenter() const
{
    //Area built cells has no points
    if(localPoints_.size()==0)
    {
        return this->center();
    }
    vector c = Zero;
    forAll(localPoints_,i)
    {
        c+=localPoints_[i];
    }

    return rotorGeometry_.cartesianToCylindrical().localPosition(c/localPoints_.size()); 
}

}

