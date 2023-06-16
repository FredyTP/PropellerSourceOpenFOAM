#include "meshCell.H"
#include "delaunayTriangulation.H"

namespace Foam
{

meshCell::meshCell(const rotorGeometry& rotorGeometry, label nBlades, const List<vector>& localVertex, const vector* center)
: gridCell(rotorGeometry), localPoints_(localVertex)
{
    if(center!=nullptr)
    {
        center_ = rotorGeometry_.cartesianToCylindrical().localPosition(*center);
    }
    else
    {
        center_ = getCellCenter();
    }

    util::geometry::sortCounterClockwise(localPoints_);
    area_ = util::geometry::poligonArea(localPoints_);
    factor_ = area_/(constant::mathematical::twoPi * center_.x());
}

vector meshCell::getCellCenter() const
{
    vector c = Zero;
    forAll(localPoints_,i)
    {
        c+=localPoints_[i];
    }

    return rotorGeometry_.cartesianToCylindrical().localPosition(c/localPoints_.size()); 
}
}

