#ifndef PROPELLER_ROTOR_DISCRETE_GRID_CELL_GRID_CELL
#define PROPELLER_ROTOR_DISCRETE_GRID_CELL_GRID_CELL

#include "scalar.H"
#include "List.H"
#include "vector.H"
#include "tensor.H"
#include "vectorField.H"
#include "scalarField.H"
#include "cylindricalCS.H"
#include "rotorGeometry.H"
#include "GeometricField.H"
#include "fvPatchField.H"
#include "volMesh.H"

namespace Foam
{
class gridCell
{
public:


    gridCell(const rotorGeometry& rotorGeometry) : rotorGeometry_(rotorGeometry), factor_(1.0) {}

    virtual ~gridCell() = default;
    /**
     * @brief Add a fvCell index and its weight to the list
     * The weight can be > 1, 
     * in the build method it is rescaled so the sum of all weights is equal 1
     * 
     * @param celli index of the fvCell
     * @param weight weight of the fvCell
     */
    void addCelli(label celli,scalar weight);

    void buildWeigths();

    void checkCells();

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
     * @brief Returns the USED CENTER in polar Coordinates z=0
     * 
     * @return const vector& 
     */
    const vector& center() const {return center_;}

    scalar radius() const {return center_.x();}

    scalar theta() const {return center_.y();}

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
    void setCenter(const vector& center);

    /**
     * @brief Get cell center in cylindrical coordinates
     * 
     * @return vector 
     */
    virtual vector getCellCenter() const = 0;


    void centerFromClosestCell(const vectorField& cellCenters);


    /**
     * FIELD RELATED FUNCTIONS *
    */

    virtual vector scaleForce(const vector& globalForce);

    
    template<class type>
    void applySource(GeometricField<type, fvPatchField, volMesh>& dstField,const scalarField& cellVol, const type& scaledSource) const;

    template<class type>
    void applySource(Field<type>& dstField, const scalarField& cellVol, const type& scaledSource) const;
    
    template<class type>
    void applyField(Field<type>& field, const type& value) const;
    
    template<class type>
    void applyField(Field<type>& field,const List<type>& values) const;
    
    label interpolatingCelli() const { return interpolatingCell_;}
private:
    vector center_; //r-theta-0
    tensor localBlade_; //(x: LE->TE) (y: radial) (z: "up")
protected:
    const rotorGeometry& rotorGeometry_;

    scalar factor_=1.0;
    List<label> cellis_;
    List<scalar> weights_;
    label interpolatingCell_=-1;

    void updateBladeTensor();
};

template <class type>
inline void gridCell::applySource(GeometricField<type, fvPatchField, volMesh> &dstField, const scalarField &cellVol, const type &scaledSource) const
{
    this->applySource<type>(dstField.ref(false),cellVol,scaledSource);
}

template <class type>
inline void gridCell::applySource(Field<type> &source, const scalarField &cellVol, const type &scaledSource) const
{
    forAll(cellis_,k)
    {
        source[cellis_[k]] =  -(weights_[k] * scaledSource/cellVol[cellis_[k]]);
    }
}
template <class type>
inline void gridCell::applyField(Field<type> &field, const type &value) const
{
    forAll(cellis_,i)
    {
        field[cellis_[i]]=value;
    }
}

template <class type>
inline void gridCell::applyField(Field<type> &field, const List<type> &values) const
{
    forAll(cellis_,i)
    {
        field[cellis_[i]]=values[i];
    }
}
}

#endif /* PROPELLER_ROTOR_DISCRETE_GRID_CELL_GRID_CELL */
