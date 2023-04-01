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
    List<point> cylPoints_;//  (Cylindrical r-theta-z)

    //- All points used in local cartesian coordinates
    List<point> carPoints_;

    List<bool> integrationPoints_;
    
    //- All cells build
    PtrList<rotorCell> rotorCells_;

    //- Total disk area
    scalar area_ = 0.0;

    //- Rotation tensor from blade local reference to global reference frame
    List<tensor> localBlade_; //(x: LE->TE) (y: radial) (z: "up")

    discreteMode discreteMode_ = discreteMode::dmNone;
public:
TypeName("rotorDiscrete");

    //Contructors
        //- Construct from dictionary
        rotorDiscrete();

        void fromRotorMesh(const rotorFvMeshSel& rotorFvMeshSel, word integration, bool correctCenters,label refinementLevel);
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

protected:
        tensor bladeLocalFromPoint(const point& localPoint) const;

};


}
#endif /* PROPELLER_ROTOR_DISCRETE_ROTOR_DISCRETE */
