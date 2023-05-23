#ifndef PROPELLER_PROPELLER_MODEL_PROPELLER_MODEL
#define PROPELLER_PROPELLER_MODEL_PROPELLER_MODEL



#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "rotorFvMeshSel.H"
#include "rotorGrid.H"
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


    void setRotorMesh(const rotorFvMeshSel* rotorFvMeshSel){ rotorFvMeshSel_ = rotorFvMeshSel;}

    virtual void build(const rotorGeometry& rotorGeometry){}

    virtual propellerResult calculate(const vectorField& U, const scalarField* rhoField, scalar angularVelocity, volVectorField& force, scalar theta){return propellerResult();}

    virtual const autoPtr<rotorGrid>& grid() = 0;

    const fvMesh& mesh() const { return rotorFvMeshSel_->mesh();}
};

}

#endif /* PROPELLER_PROPELLER_MODEL_PROPELLER_MODEL */
