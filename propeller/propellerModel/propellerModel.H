#ifndef PROPELLER_MODEL_PROPELLER_MODEL
#define PROPELLER_MODEL_PROPELLER_MODEL

#endif /* PROPELLER_MODEL_PROPELLER_MODEL */

#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "rotorMesh.H"
#include "rotorDiscrete.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class propellerModel Declaration
\*---------------------------------------------------------------------------*/

class propellerResult
{
    public:
    vector force {0,0,0}; //x: y: z: thrust in local coordinates
    scalar torque {0}; //Aerodynamic torque on rotor axis
    scalar power {0}; //Power needed to rotate the propeller
    scalar J{0}; // refV / (n*D)
    scalar eta{0}; // Thrust* refV / power

    inline scalar& thrust() {return force.z();}
    void info(word name)
    {
        Info<<name<<" result:"<<endl;
        Info<<"\tThrust = "<<thrust()<<endl;
        Info<<"\tTorque = "<<torque<<endl;
        Info<<"\tPower = "<<power<<endl;
        Info<<"\tJ = "<<J<<endl;
        Info<<"\tEta = "<<eta<<endl;
    }

};


class propellerModel
{

protected:
    const rotorMesh* rotorMesh_;
    scalar refV=1; //Refence velocity for adimensional
    scalar refRho=1; //Reference density 
public:
    TypeName("propellerModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            propellerModel,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );

    // Constructors

        //- Construct from components
        propellerModel
        (
            const dictionary& dict,
            const word& name
        );

    // Selectors

        //- Return a reference to the selected trim model
        //- automations the construction of propeller models
        //- from run time selection tables
        static autoPtr<propellerModel> New
        (
            const dictionary& dict
        );

    //Destructor

        virtual ~propellerModel() = default;

    // Abstract functions
        virtual scalar radius() const = 0;

    void setRotorMesh(const rotorMesh* rotorMesh){ rotorMesh_ = rotorMesh;}

    virtual void build(const rotorGeometry& rotorGeometry){}

    virtual propellerResult calculate(const vectorField& U, volVectorField& force){return propellerResult();}

    virtual const rotorDiscrete& rDiscrete() = 0;

    void setRefV(scalar v){refV=v;}
    void setRefRho(scalar rho){refRho = rho;}
};

}
