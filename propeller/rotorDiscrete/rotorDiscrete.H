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
       
    rotorGrid grid;

    discreteMode discreteMode_ = discreteMode::dmNone;

    rotorCell::integrationMode integrationMode_;
    discreteMethod discreteMethod_;

    bool correctCenters_ = false;
    bool includeVertex_ = true;
    label refinementLevel_;

    label nRadial;
    label nAzimutal;
    
    
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
        const rotorGeometry& geometry() const {return rotorGeometry_;}

        discreteMode mode() const {return discreteMode_;}

private:
        tensor bladeLocalFromPoint(const point& localPoint) const;

};    


}
#endif /* PROPELLER_ROTOR_DISCRETE_ROTOR_DISCRETE */
