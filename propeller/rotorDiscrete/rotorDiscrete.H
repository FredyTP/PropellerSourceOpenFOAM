#ifndef PROPELLER_ROTOR_DISCRETE_ROTOR_DISCRETE
#define PROPELLER_ROTOR_DISCRETE_ROTOR_DISCRETE

#include "fvMesh.H"
#include "cuttingPlane.H"
#include "Enum.H"
#include "cylindricalCS.H"
#include "cartesianCS.H"
#include "rotorFvMeshSel.H"
#include "rotorGeometry.H"

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

    word integrationMode_;
public:
TypeName("rotorDiscrete");

    //Contructors
        //- Construct from dictionary
        rotorDiscrete();

        /**
         * @brief Creates voronoid based mesh from cell center proyection
         * Robust and 100% tested
         * @param rotorFvMeshSel 
         * @param integration 
         * @param correctCenters 
         * @param refinementLevel 
         */
        void fromMeshVoronoid(const rotorFvMeshSel& rotorFvMeshSel, word integration, bool correctCenters,label refinementLevel);
        
        /**
         * @brief Creates rotor mesh from intersecting mesh edges with the plane
         * Not optimized, and edge intersections are computed more than once for each edge
         * There are repeated points in the list
         * 
         * Also not accounting if rotordiscrete < rotormesh, no check on centers
         * 
         * @param rotorFvMeshSel 
         * @param integration 
         * @param correctCenters 
         * @param refinementLevel 
         */
        void fromMeshIntersect(const rotorFvMeshSel &rotorFvMeshSel,word integration,bool correctCenters,label refinementLevel);
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

        word integrationMode() const { return integrationMode_;}

        void writeArea(word propName,const fvMesh& mesh) const ;
protected:
        tensor bladeLocalFromPoint(const point& localPoint) const;

};        //TODO:
        //Update:
        //Using cellweights still extreme slow because the constructor
        //Interpolate on every cell face
        //This function should return a reference to a vectorField


}
#endif /* PROPELLER_ROTOR_DISCRETE_ROTOR_DISCRETE */
