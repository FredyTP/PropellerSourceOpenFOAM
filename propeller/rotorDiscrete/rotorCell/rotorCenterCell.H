#ifndef PROPELLER_ROTOR_DISCRETE_ROTOR_CELL_ROTOR_CENTER_CELL
#define PROPELLER_ROTOR_DISCRETE_ROTOR_CELL_ROTOR_CENTER_CELL

#include "rotorCell.H"

namespace Foam
{

class rotorCenterCell : public rotorCell
{
public:
    rotorCenterCell(label center, List<label>& vertex,const List<point>& points,label cellIdx)
        : rotorCell(center,vertex,points,cellIdx), integrationPoints_({center})
    {

    }
    
    virtual vector integrateField(const List<vector>& field) const
    {
        return (field[this->center()] * this->area());
    }

    virtual const List<label>& integrationPoints() const
    {
        return integrationPoints_;
    }



private:

    List<label> integrationPoints_;

};
}


#endif /* PROPELLER_ROTOR_DISCRETE_ROTOR_CELL_ROTOR_CENTER_CELL */


