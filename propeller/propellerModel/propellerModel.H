#ifndef PROPELLER_MODEL_PROPELLER_MODEL
#define PROPELLER_MODEL_PROPELLER_MODEL

#endif /* PROPELLER_MODEL_PROPELLER_MODEL */

#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "rotorFvMeshSel.H"
#include "rotorDiscrete.H"
#include "propellerResult.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class propellerModel Declaration
\*---------------------------------------------------------------------------*/

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
