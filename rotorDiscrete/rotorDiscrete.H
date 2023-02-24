#ifndef PROPELLER_ROTOR_DISCRETE_ROTOR_DISCRETE
#define PROPELLER_ROTOR_DISCRETE_ROTOR_DISCRETE

#include "fvMesh.H"
#include "cuttingPlane.H"
#include "Enum.H"
#include "cylindricalCS.H"
#include "cartesianCS.H"
#include "rotorMesh.H"
#include "rotorGeometry.H"

namespace Foam
{

class rotorDiscrete
{

    enum discreteMode
    {
        dmNone,
        dmMesh,
        dmUniform
    };

private:

    //COORDINATE SYSTEMS
    // LabCoordinateSystem -> localCartesian -> localCylindrical
    //- Rotor local cartesian coordinate system (x-y-z)
    //Angular velocity (w->) is defined to be ( 0 0 w) in the local cartesian system, 
    //positive counter-clockwise
    coordSystem::cartesian localCS_;
    //- Rotor local cylindrical coordinate system (r-theta-z)
    //from local cartesian system
    coordSystem::cylindrical cylCS_;
    vector psiOrigin_;
    rotorGeometry rotorGeometry_;
    //- Cell centre positions in local rotor frame
    //  (Cylindrical r-theta-z)
    List<point> cylPoints_;

    //- Rotation tensor from blade local reference to global reference frame
    List<tensor> localBlade_; 

    discreteMode discreteMode_ = discreteMode::dmNone;
public:
TypeName("rotorDiscrete");

    //Contructors
        //- Construct from dictionary
        rotorDiscrete();

        void fromRotorMesh(const rotorMesh& rotorMesh);
    //Destructor
        virtual ~rotorDiscrete() = default;
    //IO
        //- Read source dictionary
        virtual bool read(const dictionary& dict);

        void buildCoordinateSystem(const rotorGeometry& geometry);

        const coordSystem::cylindrical& cylindrical() const {return cylCS_;}
        const coordSystem::cartesian& cartesian() const {return localCS_;}

        const List<point>& cylPoints() const {return cylPoints_;}

        const rotorGeometry& geometry() const {return rotorGeometry_;}
    //Getters

protected:


};


}
#endif /* PROPELLER_ROTOR_DISCRETE_ROTOR_DISCRETE */
