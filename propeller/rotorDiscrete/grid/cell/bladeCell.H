#ifndef PROPELLER_ROTOR_DISCRETE_BLADE_CELL
#define PROPELLER_ROTOR_DISCRETE_BLADE_CELL

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
class bladeCell : public gridCell
{
public:

    bladeCell(scalar radius0, scalar radius1, scalar chord0, scalar chord1);



    /**
     * @brief Reescales the weights and checks that the bladeCell contains atleast one cell
     * If the object is not used withit fv context this method doesnt need to be called
     */
    void build();


    /**
     * @brief Rotate the base points with the rotation matrix
     * 
     * @param rotation 
     */
    void setRotation(const tensor& rotation);

    /**
     * @brief Returns the actual vertex of the cell
     * 
     * @return List<vector> 
     */
    List<vector> actualPoints() const { return actualLocation_;}


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
    scalar dc() const {return dc_;}

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
    scalar chord0() const {return chord0_;}

    /**
     * @brief Second angular position
     * 
     * @return scalar 
     */
    scalar chord1() const {return chord1_;}

    /**
     * @brief Returns the tensor representing the local blade axis in global CS
     * 
     * @return const tensor& (x: LE->TE) (y: radial) (z: "up")
     */
    const tensor& localBlade() const { return localBlade_; }

    vector cartesianCenter() const;


    void centerFromClosestCell(const vectorField& cellCenters,  const coordSystem::cartesian &localCartesian);
    /**
     * @brief Set the local blade tensor in global coordinates
     *  (x: LE->TE) (y: radial) (z: "up")
     * @param localTensor 
     */
    void setLocalTensor(const tensor& localTensor);

    virtual vector scaleForce(const vector& globalForce);

    void applySource(vectorField& source,const scalarField& cellVol, vector& scaledForce) const;
    
    template<class type>
    void applyField(Field<type>& field, const type& value) const;
    
    template<class type>
    void applyField(Field<type>& field,const List<type>& values) const;
    
private:
    vector center_; //r-theta-0

    scalar radius0_; //i and i+1
    scalar radius1_; //i and i+1
    scalar chord0_; //i and i+1
    scalar chord1_; //i and i+1

    List<vector> localPoints_; //Cell base points in local cartesian coordinates, sorted couterclockwise
    List<vector> actualLocation_; //Cell position after rotation

    scalar dr_;
    scalar dc_;
};

template <class type>
inline void bladeCell::applyField(Field<type> &field, const type &value) const
{
    forAll(cellis_,i)
    {
        field[cellis_[i]]=value;
    }
}

template <class type>
inline void bladeCell::applyField(Field<type> &field, const List<type> &values) const
{
    forAll(cellis_,i)
    {
        field[cellis_[i]]=values[i];
    }
}
}

#endif /* PROPELLER_ROTOR_DISCRETE_BLADE_CELL */


