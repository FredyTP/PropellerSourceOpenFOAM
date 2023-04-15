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

/**
 * @brief Creates a mesh selection from geometry or cellSet/cellZone
 * 
 */
class rotorFvMeshSel
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
            smCylinder,
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
    
    //Selected mesh rotor geometry
    rotorGeometry meshGeometry_;
    scalar maxPointRadius=0;

    //- Zone name if imported from zone or cell set
    word cellsName_;

    bool findClosestCenter_;
    bool correctGeometry_;
    bool includeIfVertex_;

    bool built_ = false;

    //- Parallel used data
    List<label> parNcells; //num of cells in each processor
    List<List<label>> parCells; //cell index for each processor


public:
TypeName("rotorFvMeshSel");

    //Contructors
        //- Construct from parameters
        rotorFvMeshSel(const fvMesh& mesh);

    //Destructor
        virtual ~rotorFvMeshSel() = default;
    //IO
        //- Read source dictionary
        virtual bool read(const dictionary& dict);

    //Getters        
        inline const labelList& cells() const { return cells_; }
        const fvMesh& mesh() const { return mesh_;}
        inline const vector& direction() const {return meshGeometry_.direction();}
        inline const vector& center() const {return meshGeometry_.center();}
        inline scalar radius() const {return meshGeometry_.radius();}
        inline scalar maxRadius() const {return maxPointRadius;}
        inline bool isBuilt() const {return built_;}
    //BUILD SELECTION
        void build(rotorGeometry& rotorGeometry);

        void clear();

        const List<label>& coreCellNumber() const {return parNcells;}
        const List<List<label>>& coreCelli() const {return parCells;}

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

        void tryCorrectGeometry(rotorGeometry& rotorGeometry);

        void syncCellData();
};


}

#endif /* PROPELLER_ROTOR_MESH_ROTOR_FV_MESH_SEL */


