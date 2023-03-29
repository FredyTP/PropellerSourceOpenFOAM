#ifndef PROPELLER_ROTOR_MESH_ROTOR_FV_MESH_SEL
#define PROPELLER_ROTOR_MESH_ROTOR_FV_MESH_SEL
#include "fvMesh.H"
#include "cuttingPlane.H"
#include "Enum.H"
#include "cylindricalCS.H"
#include "cartesianCS.H"

#include "rotorCell.H"
#include "rotorGeometry.H"
#include "rotorCell.H"

#define NO_AREA -1
namespace Foam
{

class rotorMesh
{
public:
    //Public Enumerations

    //INFO: each enum must include a Enum class to map
    //each enum item to a string, so it can be selected from 
    //IOdicts

        //- Selection Mode
        enum selectionMode
        {
            smNone,
            smGeometry,
            smGeometry2,
            smCellSet,
            smCellZone
        };

        //- Names for selectionMode
        static const Enum<selectionMode> selectionModeNames_;
private:

    //- Reference to fluid mesh
    const fvMesh& mesh_;

    //- Mesh selection mode
    selectionMode selMode_;

    //Selected cells for rotor geometry
    //- Selected fvCells
    labelList cells_;
    //- Selected fvCells area over rotor
    scalarList area_;
    //- Sum of cells area = selected disk area
    scalar diskArea_;
    
    //- Local cartesian position
    List<vector> cellCenterLC_; //LC stands for local cartesian
    coordSystem::cartesian localCartesianCS_;

    //Rotor mesh
    List<point> rotorPoints_;
    List<rotorCell> rotorCells_;

    //Selected mesh rotor geometry
    rotorGeometry meshGeometry_;

    //- Zone name if imported from zone or cell set
    word cellsName_;

    bool findClosestCenter_;
    bool correctGeometry_;

    bool built_ = false;
public:
TypeName("rotorMesh");

    //Contructors
        //- Construct from parameters
        rotorMesh(const fvMesh& mesh);

    //Destructor
        virtual ~rotorMesh() = default;
    //IO
        //- Read source dictionary
        virtual bool read(const dictionary& dict);

    //Getters
        inline const List<rotorCell>& rotorCells() const {return rotorCells_;}
        inline const List<point>& rotorPoints() const {return rotorPoints_;}
        
        inline const labelList& cells() const { return cells_; }
        const fvMesh& mesh() const { return mesh_;}
        inline const vector& direction() const {return meshGeometry_.direction();}
        inline const vector& center() const {return meshGeometry_.center();}
        inline scalar radius() const {return meshGeometry_.radius();}
        inline scalar diskArea() const {return diskArea_;}
        inline const scalarList& areas() const {return area_;}
        inline bool isBuilt() const {return built_;}
    //BUILD SELECTION
        void build(rotorGeometry& rotorGeometry);

        void clear();

  

protected:
    //Geometry

        //- Update rotor center with closest cell center (Slow op)
        void tryUpdateCenter();
        
                
        //- Create mesh selection from center
        //and radius
        void createMeshSelection();
        void createMeshSelectionCilinder();

        //- Select cells from specified zone
        void loadMeshSelection();

        //- Compute rotor geometry from cell set
        //(find Center, find normal, find radius)
        void findRotorCenter();

        void findRotorNormal(rotorGeometry& rotorGeometry);

        void findRotorRadius();


        //- Compute cell area for selected cells
        void computeCellsArea();
        void computeCellsAreaVoronoid();

        void tryCorrectGeometry(rotorGeometry& rotorGeometry);

        void computeLocalPositions();
};


}

#endif /* PROPELLER_ROTOR_MESH_ROTOR_FV_MESH_SEL */


