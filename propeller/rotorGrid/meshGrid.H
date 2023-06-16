#ifndef PROPELLER_ROTOR_GRID_MESH_GRID
#define PROPELLER_ROTOR_GRID_MESH_GRID

#include "rotorGrid.H"
#include "Enum.H"

namespace Foam
{

class meshGrid : public rotorGrid
{
public:
    
    TypeName("meshGrid");
    enum class discreteMethod
    {
        voronoid,
        intersection,
        proyection
    };
    static const Enum<discreteMethod> discreteMethodNames_;
    meshGrid() = default;

    meshGrid
    (
        const dictionary& dict,
        const rotorGeometry& geometry,
        const rotorFvMeshSel &rotorFvMeshSel,
        const bladeModelS* bladeModel,
        scalar nBlades
    );
    void assignFvCells();
    bool read(const dictionary& dict);

    void writeArea(word propName) const;
    void writePythonPlotter();

private:
    void createMeshIntersect
    (
        List<point> &vertex, 
        List<label>& cellidx, 
        List<List<label>> &cellPoints
    );
    void createMeshVoronoid
    (
        List<point> &vertex, 
        List<point> &centers,
        List<List<label>> &cells
    );
    void createProyection();
    void fromMesh();
    
    void createFromData(List<point> &vertex, const List<point> &centers, const List<label>& cellidx, List<List<label>> &cells);    
    void selectInnerCells(List<label> &cellidx);
    bool cutWithCircle
    (
        List<point> &vertex,
        List<label> &cell,
        const std::function<bool(point)> &isInRegion,
        const std::function<void(point, point, point &, bool)> &findIntersection
    );

    void buildGrid();

    label refinementLevel_;
    bool correctCenters_;
    bool includeVertex_;

    //- Core area
    scalar area_;

    //- Total area gathering all core
    scalar totalArea_;
    discreteMethod discreteMethod_;

};

}

#endif /* PROPELLER_ROTOR_GRID_MESH_GRID */
