#ifndef PROPELLER_ROTOR_DISCRETE_ROTOR_DISCRETE
#define PROPELLER_ROTOR_DISCRETE_ROTOR_DISCRETE

#include "fvMesh.H"
#include "cuttingPlane.H"
#include "Enum.H"
#include "cylindricalCS.H"
#include "cartesianCS.H"
#include "rotorFvMeshSel.H"
#include "rotorGeometry.H"
#include "rotorCell.H"
#include "rotorGrid.H"


namespace Foam
{

class rotorDiscrete
{
public:
    enum discreteMode
    {
        dmNone, //No discrete mode yet, discretization is empty
        dmMesh, //Discretization from mesh centers
        dmPolar, //Polar discretization from radial and azimutal vectors
        dmUniform, //Uniformly distributed points
    };
    enum discreteMethod
    {
        dmeVoronoid,
        dmeIntersection,
        dmeOld
    };
    static const Enum<discreteMethod> discreteMethodNames_;

   
    rotorGrid grid;

private:

    //COORDINATE SYSTEMS
    //Two coordinate systems are used
    //- Cartesian CS from global cartesian to local cartesian (x, y, z(rotor axis))
    coordSystem::cartesian   carCS_;
    //- Cylindrical CS, from global cartesian to local cylindrical (R - theta - z)
    coordSystem::cylindrical cylCS_;
    //- Cartesian to cilindrical CS
    coordSystem::cylindrical carToCylCS_;

    //- Used rotor geometry
    rotorGeometry rotorGeometry_;
    
    //- All points used in cilindrical coordinates
    //Contains all cell vertex points(core-copied) and local core cell centers/extra points
    List<point> cylPoints_;//  (Cylindrical r-theta-z)

    //- All points used in local cartesian coordinates
    //Contains all cell vertex points(core-copied) and local core cell centers/extra points
    List<point> carPoints_;

    //- Contains which of the carPoints_ are required for integration
    List<bool> integrationPoints_;
    
    //- All cells build used in this core
    PtrList<rotorCell> rotorCells_;

    //- Core disk area
    scalar area_ = 0.0;
    scalar totalArea_ = 0.0;

    //- Rotation tensor from blade local reference to global reference frame
    List<tensor> localBlade_; //(x: LE->TE) (y: radial) (z: "up")

    discreteMode discreteMode_ = discreteMode::dmNone;

    rotorCell::integrationMode integrationMode_;
    discreteMethod discreteMethod_;

    bool correctCenters_ = false;
    bool includeVertex_ = true;
    label refinementLevel_;
    
    
public:
TypeName("rotorDiscrete");

    //Contructors
        //- Construct from dictionary
        rotorDiscrete(const dictionary& dict);


    //Builders
        /**
         * @brief Construct rotor discre from mesh
         * 
         * @param rotorFvMeshSel 
         */
        void fromMesh
        (
            const rotorFvMeshSel &rotorFvMeshSel
        );

    //Output
        /**
         * @brief creates a python file to plot the discrete
         * 
         * @param process 
         */
        void writePythonPlotter(word process="");
    //Destructor
        virtual ~rotorDiscrete() = default;
    //IO
        //- Read source dictionary
        virtual bool read(const dictionary& dict);

        void buildCoordinateSystem(const rotorGeometry& geometry);

     //Getters
        const coordSystem::cylindrical& cylindrical() const {return cylCS_;}
        const coordSystem::cartesian& cartesian() const {return carCS_;}
        const List<point>& cylPoints() const {return cylPoints_;}

        const rotorGeometry& geometry() const {return rotorGeometry_;}
        
        const List<tensor>& localBladeCS() const {return localBlade_;}

        label size() const {return cylPoints_.size();}

        const PtrList<rotorCell>& rotorCells() const {return rotorCells_;}

        discreteMode mode() const {return discreteMode_;}

        const List<bool>& integrationPoints() const { return integrationPoints_; }

        rotorCell::integrationMode integrationMode() const { return integrationMode_;}

        void writeArea(word propName,const fvMesh& mesh) const ;

private:

        void createMeshOld
        (
            const rotorFvMeshSel& rotorFvMeshSel
        );
        /**
         * @brief Create a Mesh Voronoid object, 
         * needs to gather all points from all the cores
         * 
         * @param rotorFvMeshSel 
         * @param vertex 
         * @param centers 
         * @param cells 
         */
        void createMeshVoronoid
        (
            List<point> &vertex, 
            List<point> &centers,
            List<List<label>> &cells
        );
        
        /**
         * @brief Create a Mesh Intersect object
         * Obtain a list of cell points indexing to carPoints_
         * 
         * @param rotorFvMeshSel 
         * @param cellPoints 
         */
        void createMeshIntersect
        (
            const rotorFvMeshSel& rotorFvMeshSel, 
            List<point> &vertex, 
            List<label>& cellidx, 
            List<List<label>> &cellPoints
        );

        /**
         * @brief Select cells inside the rotor (center or vertex based)
         * 
         * @param rotorFvMeshSel 
         * @param cellidx 
         */
        void selectInnerCells
        (
            const rotorFvMeshSel& rotorFvMeshSel, 
            List<label>& cellidx
        );

        /**
         * @brief 
         * 
         * @param vertex 
         * @param cell 
         * @param isInRegion 
         * @param findIntersection 
         * @return true 
         * @return false 
         */
        bool cutWithCircle
        (
            List<point> &vertex,
            List<label> &cell,
            const std::function<bool(point)> &isInRegion,
            const std::function<void(point, point, point &, bool)> &findIntersection
        );
        /**
         * @brief Create a From Data object
         * 
         * @param vertex 
         * @param centers 
         * @param cellidx 
         * @param cells 
         */
        void createFromData
        (
            const List<point> &vertex, 
            const List<point> &centers,
            const List<label>& cellidx, 
            List<List<label>> &cells
        );
       
public:
        tensor bladeLocalFromPoint(const point& localPoint) const;

};        //TODO:
        //Update:
        //Using cellweights still extreme slow because the constructor
        //Interpolate on every cell face
        //This function should return a reference to a vectorField


}
#endif /* PROPELLER_ROTOR_DISCRETE_ROTOR_DISCRETE */
