#ifndef PROPELLER_ROTOR_GRID_CELL_MESH_CELL
#define PROPELLER_ROTOR_GRID_CELL_MESH_CELL

#include "scalar.H"
#include "List.H"
#include "vector.H"
#include "tensor.H"
#include "vectorField.H"
#include "scalarField.H"
#include "cylindricalCS.H"

#include "cartesianCS.H"
#include "gridCell.H"

namespace Foam
{
class meshCell : public gridCell
{
public:
    meshCell(const rotorGeometry& rotorGeometry, label nBlades, const List<vector>& localVertex, const vector* center);


    /**
     * @brief Returns the radial size of the cell (radius[i+1]-radius[i])
     * 
     * @return scalar 
     */
    scalar dr() const { return dr_;}

    /**
     * @brief Get cell center in cylindrical coordinates
     * 
     * @return vector 
     */
    vector getCellCenter() const;

private:
    vector center_; //r-theta-0

    List<vector> localPoints_; //Cell base points in local cartesian coordinates, sorted couterclockwise

    scalar dr_;
    scalar dc_;
    scalar area_;
};


}


#endif /* PROPELLER_ROTOR_GRID_CELL_MESH_CELL */
