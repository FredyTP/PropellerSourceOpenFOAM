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

/**
 * Propeller result integrated data information
 * Adimensional parameters defined as website : 
 * https://m-selig.ae.illinois.edu/props/propDB.html
 * suggests
*/
class propellerResult
{
    public:
    vector force {0,0,0}; //- x: y: z: thrust in local coordinates
    scalar torque {0};    //- Aerodynamic torque on rotor axis
    scalar power {0};     //- Power needed to rotate the propeller
    scalar J{0};          //- Adimensional: refV / (n*D)
    scalar eta{0};        //- Thrust* refV / power
    scalar CT{0};         //- Thrust coefficient T/(rho*n^2*D^4)
    scalar CP{0};         //- Power coefficient Power/(rho*n^3*D^5)

    inline scalar thrust() const {return force.z();}
    void info(word name)
    {
        Info<<name<<" result:"<<endl;
        Info<<"\tThrust = "<<thrust()<<endl;
        Info<<"\tTorque = "<<torque<<endl;
        Info<<"\tPower = "<<power<<endl;
        Info<<"\tJ = "<<J<<endl;
        Info<<"\tEta = "<<eta<<endl;
    }

    /**
     * Update Js with force and power already set
     * J defined as Adimensional: refV / (n*D)
    */
    void updateJ(scalar vRef, scalar omega, scalar radius)
    {
        scalar vRot = omega/constant::mathematical::pi * radius;
        if(vRot == 0)
        {
            J = 0;
        }
        else
        {
            J = vRef/vRot;
        }
    }

    /**
     * @brief: Update eta value with force and power already set
     * eta defined as Adimensional: Thrust * vRef / power
    */
    void updateEta(scalar vRef)
    {
        if(power == 0)
        {
            eta = 0;
        }
        else
        {
            eta=thrust() * vRef / power;
        }
    }

    void updateCT(scalar rhoRef,scalar omega, scalar radius)
    {
        scalar den = rhoRef * 
                pow((omega / constant::mathematical::piByTwo),2)*
                pow(2*radius,4);
        if(den == 0)
        {
            CT=0;
        }
        else
        {
            CT = thrust()/den;
        }
    }

    void updateCP(scalar rhoRef, scalar omega,scalar radius)
    {
        scalar den = rhoRef * 
                pow((omega / constant::mathematical::piByTwo),3)*
                pow(2*radius,5);
        if(den == 0)
        {
            CP=0;
        }
        else
        {
            CP = power/den;
        }
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

    virtual propellerResult calculate(const vectorField& U,scalar angularVelocity, volVectorField& force){return propellerResult();}

    virtual const rotorDiscrete& rDiscrete() = 0;

    void setRefV(scalar v){refV=v;}
    void setRefRho(scalar rho){refRho = rho;}
};

}
