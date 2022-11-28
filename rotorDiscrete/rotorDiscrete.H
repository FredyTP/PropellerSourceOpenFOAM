#ifndef PROPELLER_ROTOR_GEOMETRY_ROTOR_GEOMETRY
#define PROPELLER_ROTOR_GEOMETRY_ROTOR_GEOMETRY

#include "fvMesh.H"
#include "cuttingPlane.H"
#include "Enum.H"
#include "cylindricalCS.H"

namespace Foam
{

class rotorDiscrete
{
public:
    //Public Enumerations

    //INFO: each enum must include a Enum class to map
    //each enum item to a string, so it can be selected from 
    //IOdicts

        //- Selection Mode
        enum selectionMode
        {
            smGeometry,
            smCellSet
        };

        //- Names for selectionMode
        static const Enum<selectionMode> selectionModeNames_;
private:
    const fvMesh& mesh_;

    selectionMode selMode_;

    //Selected cells for rotor geometry
    labelList cells_;
    scalarList area_;
    scalar diskArea_;

    //COORDINATE SYSTEMS
    //- Rotor local cylindrical coordinate system (r-theta-z)
    coordSystem::cylindrical coordSys_;
    vector psiOrigin_;
    //- Cell centre positions in local rotor frame
    //  (Cylindrical r-theta-z)
    List<point> cylCellCenter;

    //Specified atributes
    scalar radius_;
    vector rotorCenter_;
    vector rotorDir_;
    bool isClosestCenter_;

    //Parameters obtained after cell selection
    vector realCenter_;
public:
TypeName("rotorDiscrete");

    //Contructors
        //- Construct from dictionary
        rotorDiscrete(scalar radius, const fvMesh& mesh, const dictionary& dict);

    //Destructor
        virtual ~rotorDiscrete() = default;
    //IO
        //- Read source dictionary
        virtual bool read(const dictionary& dict);

    //Getters

        inline const labelList& cells(){ return cells_; }
        inline const vector& direction(){return rotorDir_;}
        inline const vector& center(){return realCenter_;}

protected:
    //Geometry

        void updateCenter();
        
        //- Create mesh selection from center
        //and radius
        void createMeshSelection();

        void findArea();



};


}
#endif /* PROPELLER_ROTOR_GEOMETRY_ROTOR_GEOMETRY */
