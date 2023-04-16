#ifndef PROPELLER_ROTOR_DISCRETE_ROTOR_CELL_ROTOR_CELL
#define PROPELLER_ROTOR_DISCRETE_ROTOR_CELL_ROTOR_CELL

#include "point.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "Enum.H"
namespace Foam
{


class rotorCell
{
public:

    enum integrationMode
    {
        imCenter,
        imTri,
        imTriCenter
    };
    static const Enum<integrationMode> integrationModeNames_;

    rotorCell(label center, List<label>& vertex,const List<point>& points,label cellIdx)
        :center_(center), vertex_(vertex), pointData_(&points), celli_(cellIdx)
    {
        this->computeArea();
    }
    
    virtual ~rotorCell() = default;

    // -----PURE VIRTUAL METHODS  
    

    /**
     * @brief: Integrates field over the cells integration points.
     * 
     * @details:
     * The points obtanined from function "integrationPoints()" are
     * guaranteed included in the field in the correct index
     * 
     * Accessing field[integrationPoint[0]] -> returns the field value at integrationPoint 0
    */
    virtual vector integrateField(const List<vector>& field) const = 0;


    /**
     * @brief: Returns a list of integration points index
    */
    virtual const List<label>& integrationPoints() const = 0;


    //-------------------------------------------------------------//

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

    inline const point& position(label i) const { return (*pointData_)[i];} 
    inline const point& vertexPosition(label vertexi) const { return position(vertex_[vertexi]);}
    inline const point& centerPosition() const { return position(center_);}
    inline label nVertex() const {return vertex_.size();}
    inline const List<label>& vertex() const {return vertex_;}
    inline label center() const { return center_;}
             
    static autoPtr<rotorCell> New(integrationMode type, label center, List<label>& vertex, const List<point> &points, label cellIdx);

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
            area_ += (vertexPosition(j).x() + vertexPosition(i).x()) * (vertexPosition(j).y() - vertexPosition(i).y());
            j = i;
        }
        area_ = std::abs(area_/ 2.0);
    }


private:
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


