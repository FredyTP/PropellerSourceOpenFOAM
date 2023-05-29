#ifndef PROPELLER_PROPELLER_MODEL_PROPELLER_MODEL
#define PROPELLER_PROPELLER_MODEL_PROPELLER_MODEL



#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "rotorFvMeshSel.H"
#include "rotorGrid.H"
#include "propellerResult.H"
#include "rotorControl.H"
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
    virtual void nextTimeStep(scalar dt){}
    virtual propellerResult calculate(const vectorField& U, const scalarField* rhoField, volVectorField& force){return propellerResult();}
    virtual propellerResult calculate(const vectorField& U, const scalarField* rhoField) const {return propellerResult();}

    virtual const autoPtr<rotorGrid>& grid() = 0;
    virtual rotorControl& control() = 0;

    const fvMesh& mesh() const { return rotorFvMeshSel_->mesh();}
    
    static tensor bladeTensor(const coordSystem::cylindrical &cylCS, const point &localPoint, scalar flapping, scalar sweep);

};

}

#endif /* PROPELLER_PROPELLER_MODEL_PROPELLER_MODEL */
