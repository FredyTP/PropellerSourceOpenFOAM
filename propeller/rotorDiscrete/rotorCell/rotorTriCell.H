#ifndef PROPELLER_ROTOR_DISCRETE_ROTOR_CELL_ROTOR_TRI_CELL
#define PROPELLER_ROTOR_DISCRETE_ROTOR_CELL_ROTOR_TRI_CELL

#include "List.H"
#include "vector.H"
#include "delaunayTriangulation.H"
#include "rotorCell.H"

namespace Foam
{
    /**
     * Cell for calculation and integration of pressure
     * Use diferent list for integration points rather than using the same as vertex
    */
class rotorTriCell : public rotorCell
{
public:

    rotorTriCell(label center, List<label>& vertex,const List<point>& points,label cellIdx, bool useCenter)
        : rotorCell(center,vertex,points,cellIdx)
    {

        integrationPoints_ = this->vertex();
        if(useCenter)
        {
            integrationPoints_.append(this->center());
        }

        List<point> pts(integrationPoints_.size());
        forAll(pts,i)
        {
            pts[i] = this->position(integrationPoints_[i]);
        }


        delaunayTriangulation::Triangulate(pts,tri_);

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
        triTotalArea_ = areaTri;
    }
    

    virtual vector integrateField(const List<vector>& field) const
    {
        vector value{0,0,0};

        forAll(tri_,i)
        {
            value += (triArea_[i] * (field[tri_[i][0]]+field[tri_[i][1]]+field[tri_[i][2]])/3.0);
        }
        
        return value;
    }

    virtual const List<label>& integrationPoints() const
    {
        return integrationPoints_;
    }
  

    //GETTERS
    inline const List<FixedList<label,3>>& tri() const { return tri_;}
    inline const List<scalar>& triArea() const { return triArea_;}
    inline scalar areaTri() const { return triTotalArea_; }
    

private:
    label globalIndex(label localIndex) const {return integrationPoints_[localIndex];}
    
    /**
     * Comput tri area
    */
     scalar computeAreaTri(const FixedList<label,3>& t)
    {
        scalar area = 0.0;

        label j = t.size() - 1;
        for (label i = 0; i < t.size(); i++)
        {
            area += (this->position(t[j]).x() + this->position(t[i]).x()) * (this->position(t[j]).y() - this->position(t[i]).y());
            j = i;
        }
        return std::abs(area / 2.0);
    }

    List<label> integrationPoints_;
    
    List<FixedList<label,3>> tri_;
    List<scalar> triArea_;
    scalar triTotalArea_;
};
}
#endif /* PROPELLER_ROTOR_DISCRETE_ROTOR_CELL_ROTOR_TRI_CELL */

