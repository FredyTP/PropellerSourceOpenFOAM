#ifndef PROPELLER_ROTOR_DISCRETE_ROTOR_CELL_ROTOR_CELL
#define PROPELLER_ROTOR_DISCRETE_ROTOR_CELL_ROTOR_CELL

#include "dictionary.H"
#include "autoPtr.H"

namespace Foam
{


class rotorCell
{
public:

    rotorCell(label center, List<label> vertex,const List<point>& points,label cellIdx)
        :center_(center), vertex_(vertex), pointData_(&points), celli_(cellIdx)
    {
        this->computeArea();
    }
    

    /**
     * @brief: Integrates field over the cells integration points.
     * 
     * @details:
     * The points obtanined from function "integrationPoints()" are
     * guaranteed included in the field in the correct index
     * 
     * Accessing field[integrationPoint[0]] -> returns the field value at integrationPoint 0
    */
    template<typename fieldType>
    virtual scalar integrateField(const List<fieldType>& field) const = 0;


    /**
     * @brief: Returns a list of integration points index
    */
    virtual const List<label>& integrationPoints() const = 0;

    /**
     * @brief: Updates the integration points values needed
    */
    void updateIntegrationList(List<bool>& intList) const
    {
        const auto& intPoi = this->integrationPoints();

        forAll(intPoi,i)
        {
            intList[intPoi[i]] = true;
        }
    }

    // GETTERS
    inline scalar area() const { return area_; }
    inline scalar celli() const {return celli_; }

    inline const point& position(label i) const { return (*pointData)[i];} 
    inline const point& vertexPosition(label vertexi) const { return postion(vertex_[i]);}

protected:

    /**
     * Compute cell area
    */
    void computeArea()
    {
        area_ = 0.0;

        label j = vertex_.size() - 1;
        for (label i = 0; i < vertex_.size(); i++)
        {
            area += (position(vertex_[j]).x() + position(vertex_[i]).x()) * (position(vertex_[j]).y() - position(vertex_[i]).y());
            j = i;
        }
        std::abs(area_/ 2.0);
    }


    //Idx of center position   
    label center_; 

    //List of vertex idx sortex counterclockwise
    List<label> vertex_;

    //Points array ptr
    const List<point>* pointData_;
    
    //Cell idx of the fvMesh cell
    label celli_;
    
    //Area of the cell
    scalar area_;

};



}

#endif /* PROPELLER_ROTOR_DISCRETE_ROTOR_CELL_ROTOR_CELL */
