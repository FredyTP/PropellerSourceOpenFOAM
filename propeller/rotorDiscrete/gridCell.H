#ifndef PROPELLER_ROTOR_DISCRETE_GRID_CELL
#define PROPELLER_ROTOR_DISCRETE_GRID_CELL

#include "scalar.H"
#include "List.H"
#include "vector.H"

namespace Foam
{
class gridCell
{
public:
    /**
     * @brief Construct a new grid Cell object
     * 
     * @param radius0 radius where cell begins
     * @param radius1 radius where cell ends
     * @param theta0 theta where cell begins [-pi, pi]
     * @param theta1 theta where cell ends [-pi, pi]
     */
    gridCell(scalar radius0,scalar radius1, scalar theta0, scalar theta1);


    /**
     * @brief Add a fvCell index and its weight to the list
     * The weight can be > 1, 
     * in the build method it is rescaled so the sum of all weights is equal 1
     * 
     * @param celli index of the fvCell
     * @param weight weight of the fvCell
     */
    void addCelli(label celli,scalar weight);

    /**
     * @brief Reescales the weights and checks that the gridCell contains atleast one cell
     * If the object is not used withit fv context this method doesnt need to be called
     */
    void build();

    /**
     * @brief Returns a list of the fvCell indexs
     * 
     * @return const labelList& 
     */
    const labelList& cellis() const { return cellis_; }

    /**
     * @brief Returns a list of the fvCell weights
     * 
     * @return const labelList& 
     */
    const List<scalar>& weights() const {return weights_;}

    /**
     * @brief Returns the center radius
     * 
     * @return scalar 
     */
    scalar radius() const {return center_.x();}

    /**
     * @brief Returns the center theta
     * 
     * @return scalar 
     */
    scalar theta()  const {return center_.y();}

    /**
     * @brief Returns the center in polar Coordinates z=0
     * 
     * @return const vector& 
     */
    const vector& center() const {return center_;}

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

    /**
     * @brief Returns the tensor representing the local blade axis in global CS
     * 
     * @return const tensor& (x: LE->TE) (y: radial) (z: "up")
     */
    const tensor& localBlade() const { return localBlade_; }

    /**
     * @brief Set the cell center in local polar coordinates z=0
     * 
     * @param center 
     */
    void setCenter(vector& center);


    void centerFromClosestCell(const vectorField& cellCenters, const coordSystem::cylindrical& localCyl);
    /**
     * @brief Set the local blade tensor in global coordinates
     *  (x: LE->TE) (y: radial) (z: "up")
     * @param localTensor 
     */
    void setLocalTensor(tensor& localTensor);


    label interpolatingCelli() const { return interpolatingCell_;}
private:
    vector center_; //r-theta-0

    scalar radius0_; //i and i+1
    scalar radius1_; //i and i+1
    scalar theta0_; //i and i+1
    scalar theta1_; //i and i+1
    
    tensor localBlade_; //(x: LE->TE) (y: radial) (z: "up")

    scalar dr_;
    scalar dt_;

    List<label> cellis_;
    List<scalar> weights_;

    label interpolatingCell_=-1;
};

}

#endif /* PROPELLER_ROTOR_DISCRETE_GRID_CELL */
