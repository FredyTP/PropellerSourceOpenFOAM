#ifndef PROPELLER_MODEL_PROPELLER_MODEL
#define PROPELLER_MODEL_PROPELLER_MODEL

#endif /* PROPELLER_MODEL_PROPELLER_MODEL */

#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "rotorFvMeshSel.H"
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
    vector torque {0,0,0};    //- Aerodynamic torque on rotor axis
    scalar power {0};     //- Power needed to rotate the propeller
    scalar J{0};          //- Adimensional: refV / (n*D)
    scalar eta{0};        //- Thrust* refV / power
    scalar CT{0};         //- Thrust coefficient T/(rho*n^2*D^4)
    scalar CP{0};         //- Power coefficient Power/(rho*n^3*D^5)

    propellerResult()
    {
        propellerResult::OutputDefinition(Info)<<endl;
    }
    inline scalar thrust() const {return force.z();}

    static Ostream& OutputDefinition(Ostream& os)
    {
        if(!definitionShown_)
        {
            definitionShown_ = true;
            indent(os) << "Rotor adimiensional parameters definitions: "<<endl;
            os.incrIndent();
            indent(os) <<"- J (Advance Coefficient):  Vref/(n*D)"<<endl;
            indent(os) <<"- Eta (Efficiency): (T * Vref)/P"<<endl;
            indent(os) <<"- CT (Thrust Coefficient): T/(rho * n^2 * D^4)"<<endl;
            indent(os) <<"- CP (Power Coefficient): P/(rho * n^3 * D^5)"<<endl;
            os.decrIndent();
            indent(os) << "Where: "<<endl;
            os.incrIndent();
            indent(os)<< "- Vref is reference velocity"<<endl;
            indent(os)<< "- n is rotational speed in revolution/second"<<endl;
            indent(os)<< "- D is rotor diameter"<<endl;
            indent(os)<< "- T is rotor thrust"<<endl;
            indent(os)<< "- P is rotor power"<<endl;
            indent(os)<< "- rho is reference density"<<endl;
            os.decrIndent();
        }
        return os;
    }

    friend Ostream& operator<<(Ostream& os, const propellerResult& pr)
    {   
        os.incrIndent();
        indent(os) <<"- Force = "<<pr.force<<endl;
        indent(os) <<"- Torque = "<<pr.torque<<endl;
        indent(os) <<"- Power = "<<pr.power<<endl;
        indent(os) <<"- J = "<<pr.J<<endl;
        indent(os) <<"- Eta = "<<pr.eta<<endl;
        indent(os) <<"- CT = "<<pr.CT<<endl;
        indent(os) <<"- CP = "<<pr.CP<<endl; 
        os.decrIndent();
        return os;
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

private:
    static bool definitionShown_;
};




class propellerModel
{

protected:
    const rotorFvMeshSel* rotorFvMeshSel_;
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

    void setRotorMesh(const rotorFvMeshSel* rotorFvMeshSel){ rotorFvMeshSel_ = rotorFvMeshSel;}

    virtual void build(const rotorGeometry& rotorGeometry){}

    virtual propellerResult calculate(const vectorField& U,scalar angularVelocity, volVectorField& force){return propellerResult();}

    virtual const rotorDiscrete& rDiscrete() = 0;

    void setRefV(scalar v){refV=v;}
    void setRefRho(scalar rho){refRho = rho;}
};

}
