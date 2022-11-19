#ifndef PROPELLER_ROTOR_GEOMETRY_ROTOR_GEOMETRY
#define PROPELLER_ROTOR_GEOMETRY_ROTOR_GEOMETRY

#include "fvMesh.H"
#include "cuttingPlane.H"

namespace Foam
{

class rotorGeometry
{
public:
    TypeName("rotorGeometry");

    //Contructors
        //- Construct from dictionary
        rotorGeometry(const fvMesh& mesh, const dictionary& dict);

    //Destructor
        virtual ~rotorGeometry() = default;
    //IO
        //- Read source dictionary
        virtual bool read(const dictionary& dict);

    //Getters

        inline const labelList& cells(){ return cells_; }
        inline const vector& direction(){return rotorDir_;}
        inline const vector& center(){return realCenter_;}
private:
    const fvMesh& mesh_;

    //Selected cells for rotor geometry
    labelList cells_;
    
    //Specified atributes
    scalar radius_;
    vector rotorCenter_;
    vector rotorDir_;
    bool isClosestCenter_;

    //Parameters obtained after cell selection
    vector realCenter_;

protected:
    //Geometry

        void updateCenter();
        
        //- Create mesh selection from center
        //and radius
        void createMeshSelection();



};


}
#endif /* PROPELLER_ROTOR_GEOMETRY_ROTOR_GEOMETRY */
