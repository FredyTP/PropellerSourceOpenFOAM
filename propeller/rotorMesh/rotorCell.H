#ifndef PROPELLER_ROTOR_MESH_ROTOR_CELL
#define PROPELLER_ROTOR_MESH_ROTOR_CELL

#include "List.H"
#include "vector.H"
#include "delaunayTriangulation.H"

namespace Foam
{
class rotorCell
{
public:

    rotorCell()
    {
        
    }
    rotorCell(label center, List<label> vertex,const List<point>& points,scalar vol)
        : center_(center), vertex_(vertex), pointData_(&points), volume_(vol)
    {
        List<point> pts(vertex_.size()+1);
        forAll(vertex_,i)
        {
            pts[i] = localPosition(i);
        }
        pts[vertex_.size()]=localPosition(vertex_.size());
        delaunayTriangulation::Triangulate(pts,tri_);
        area_ = computeArea(vertex_);
        triArea_.resize(tri_.size());

        forAll(tri_,i)
        {
            triArea_[i] = computeAreaTri(tri_[i]);
        }
    }
    const point& localPosition(label i) const {return i==vertex_.size()? (*pointData_)[center_] : (*pointData_)[vertex_[i]];}
    const point& position(label i) const {return (*pointData_)[i];}

    scalar area() { return area_;}


    //GETTERS
    inline const List<FixedList<label,3>>& tri() const { return tri_;}
    inline const List<scalar>& triArea() const { return triArea_;}
    inline label center() const { return center_;}
    inline const List<label>& vertex() const {return vertex_;}
    inline scalar volume() const {return volume_;}

private:

    /**
     * Compute cell area
    */
    scalar computeArea(const List<label>& v)
    {
        scalar area = 0.0;

        label j = v.size() - 1;
        for (label i = 0; i < v.size(); i++)
        {
            area += (localPosition(j).x() + localPosition(i).x()) * (localPosition(j).y() - localPosition(i).y());
            j = i;
        }
        return std::abs(area/ 2.0);
    }

    /**
     * Comput tri area
    */
     scalar computeAreaTri(const FixedList<label,3>& t)
    {
        scalar area = 0.0;

        label j = t.size() - 1;
        for (label i = 0; i < t.size(); i++)
        {
            area += (localPosition(t[j]).x() + localPosition(t[i]).x()) * (localPosition(t[j]).y() - localPosition(t[i]).y());
            j = i;
        }
        return std::abs(area/ 2.0);
    }

    
    List<label> vertex_;
    label center_; 
    scalar area_;
    scalar volume_;
    const List<point>* pointData_;
    List<FixedList<label,3>> tri_;
    List<scalar> triArea_;
};
}

#endif /* PROPELLER_ROTOR_MESH_ROTOR_CELL */
