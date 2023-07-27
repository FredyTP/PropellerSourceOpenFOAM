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
    /**
     * @brief Construct a new mesh Cell object from list of points
     * 
     * @param rotorGeometry 
     * @param nBlades 
     * @param localVertex 
     * @param center 
     */
    meshCell(const rotorGeometry& rotorGeometry, label celli, label nBlades, const List<vector>& localVertex, const vector* center);

    /**
     * @brief Construct a new mesh Cell object from indexed points
     * 
     * @param rotorGeometry 
     * @param nBlades 
     * @param points 
     * @param vertexIndex 
     * @param center 
     */
    meshCell(const rotorGeometry& rotorGeometry, label celli, label nBlades, const List<point>& points, const List<label>& vertexIndex, const vector* center);

    /**
     * @brief Construct a new mesh Cell object from area and center
     * 
     * @param rotorGeometry 
     * @param nBlades 
     * @param area 
     * @param center 
     */
    meshCell(const rotorGeometry& rotorGeometry, label celli, label nBlades, scalar area, const vector& center);

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

    void setCell(label celli);

    const List<vector>& localCartesianPoints() const { return localPoints_; }
private:
    List<vector> localPoints_; //Cell base points in local cartesian coordinates, sorted couterclockwise

    scalar dr_;
};


}


#endif /* PROPELLER_ROTOR_GRID_CELL_MESH_CELL */
