#ifndef PROPELLER_ROTOR_DISCRETE_ROTOR_CELL_ROTOR_TRI_CELL
#define PROPELLER_ROTOR_DISCRETE_ROTOR_CELL_ROTOR_TRI_CELL

#include "List.H"
#include "vector.H"
#include "delaunayTriangulation.H"

namespace Foam
{
    /**
     * Cell for calculation and integration of pressure
     * Use diferent list for integration points rather than using the same as vertex
    */
class rotorTriCell
{
public:

    rotorTriCell()
    {
        
    }
    rotorTriCell(label center, List<label> vertex,const List<point>& points,label cellIdx)
        : center_(center), vertex_(vertex), pointData_(&points), celli_(cellIdx)
    {
        vertex_.append(center_);
        List<point> pts(vertex_.size());
        forAll(vertex_,i)
        {
            pts[i] = localPosition(i);
        }
        //pts[vertex_.size()]=localPosition(vertex_.size());

        delaunayTriangulation::Triangulate(pts,tri_);
        area_ = computeArea(vertex_);
        triArea_.resize(tri_.size());

        scalar areaTri = 0;
        forAll(tri_,i)
        {
            tri_[i][0] = globalIndex(tri_[i][0]);
            tri_[i][1] = globalIndex(tri_[i][1]);
            tri_[i][2] = globalIndex(tri_[i][2]);

            triArea_[i] = computeAreaTri(tri_[i]);
            areaTri += triArea_[i];
        }
    }
    const label globalIndex(label localIndex) {return localIndex==vertex_.size()? center_ : vertex_[localIndex];}
    const point& localPosition(label i) const {return i==vertex_.size()? (*pointData_)[center_] : (*pointData_)[vertex_[i]];}
    const point& position(label i) const {return (*pointData_)[i];}

    scalar area() const { return area_;}

    //GETTERS
    inline const List<FixedList<label,3>>& tri() const { return tri_;}
    inline const List<scalar>& triArea() const { return triArea_;}
    inline label center() const { return center_;}
    inline const List<label>& vertex() const {return vertex_;}
    inline label celli() const {return celli_;}

private:

    
    /**
     * Comput tri area
    */
     scalar computeAreaTri(const FixedList<label,3>& t)
    {
        scalar area = 0.0;

        label j = t.size() - 1;
        for (label i = 0; i < t.size(); i++)
        {
            area += (position(t[j]).x() + position(t[i]).x()) * (position(t[j]).y() - position(t[i]).y());
            j = i;
        }
        return std::abs(area/ 2.0);
    }

    List<label> vertex_;
    label center_; 
    scalar area_;
    label celli_;
    const List<point>* pointData_;
    List<FixedList<label,3>> tri_;
    List<scalar> triArea_;
};
}
#endif /* PROPELLER_ROTOR_DISCRETE_ROTOR_CELL_ROTOR_TRI_CELL */
