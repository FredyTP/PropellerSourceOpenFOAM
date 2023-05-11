#ifndef PROPELLER_ROTOR_DISCRETE_GRID_CELL_POLAR_CELL
#define PROPELLER_ROTOR_DISCRETE_GRID_CELL_POLAR_CELL

#include "scalar.H"
#include "List.H"
#include "vector.H"
#include "tensor.H"
#include "vectorField.H"
#include "scalarField.H"
#include "cylindricalCS.H"
#include "gridCell.H"

namespace Foam
{
class polarCell : public gridCell
{
public:

    /**
     * @brief Construct a new grid Cell object
     * @param rotorGeometry geometry rotor (for CS)
     * @param radius0 radius where cell begins
     * @param radius1 radius where cell ends
     * @param theta0 theta where cell begins [-pi, pi]
     * @param theta1 theta where cell ends [-pi, pi]
     */
    polarCell(const rotorGeometry& rotorGeometry,scalar radius0,scalar radius1, scalar theta0, scalar theta1,label nBlades);



    /**
     * @brief Get cell center in cylindrical coordinates
     * 
     * @return vector 
     */
    vector getCellCenter() const;

    /**
     * @brief Returns the radial size of the cell (radius[i+1]-radius[i])
     * 
     * @return scalar 
     */
    scalar dr() const { return dr_;}

    /**
     * @brief Returns the angular size of the cell (theta[i+1]-theta[i])
     * 
     * @return scalar 
     */
    scalar dt() const {return dt_;}

    /**
     * @brief First radial position
     * 
     * @return scalar 
     */
    scalar radius0() const {return radius0_;}

    /**
     * @brief Second radial position
     * 
     * @return scalar 
     */
    scalar radius1() const {return radius1_;}

    /**
     * @brief First angular position
     * 
     * @return scalar 
     */
    scalar theta0() const {return theta0_;}

    /**
     * @brief Second angular position
     * 
     * @return scalar 
     */
    scalar theta1() const {return theta1_;}



private:

    scalar radius0_; //i and i+1
    scalar radius1_; //i and i+1
    scalar theta0_; //i and i+1
    scalar theta1_; //i and i+1
    
    scalar dr_;
    scalar dt_;



};

}

#endif /* PROPELLER_ROTOR_DISCRETE_GRID_CELL_POLAR_CELL */


