#ifndef PROPELLER_PROPELLER_MODEL_PROPELLER_MODEL
#define PROPELLER_PROPELLER_MODEL_PROPELLER_MODEL



#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "rotorFvMeshSel.H"
#include "rotorGrid.H"
#include "propellerResult.H"
#include "rotorControl.H"
#include "outputFields.H"

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
                word sourceName,
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (sourceName,dict,mesh)
        );

    // Constructors

        //- Construct from components
        propellerModel
        (
            word sourceName,
            const dictionary& dict,
            const fvMesh& mesh
        );

    
    // Selectors

        //- Return a reference to the selected trim model
        //- automations the construction of propeller models
        //- from run time selection tables
        static autoPtr<propellerModel> New
        (
            word sourceName,
            const dictionary& dict,
            const fvMesh& mesh
        );

    //Destructor

        virtual ~propellerModel() = default;


    void setRotorMesh(const rotorFvMeshSel* rotorFvMeshSel){ rotorFvMeshSel_ = rotorFvMeshSel;}

    virtual void build(const rotorGeometry& rotorGeometry){}
    virtual bool nextTimeStep(scalar dt){return false;}
    virtual propellerResult calculate(const vectorField& U, const scalarField* rhoField, volVectorField& force){return propellerResult();}
    virtual propellerResult calculate(const vectorField& U, const scalarField* rhoField) const {return propellerResult();}

    virtual const autoPtr<rotorGrid>& grid() const = 0;
    
    const fvMesh& mesh() const { return mesh_;}
    const word& sourceName() const {return sourceName_;}
    
    static tensor bladeTensor(const coordSystem::cylindrical &cylCS, const point &localPoint, scalar flapping, scalar sweep);

    virtual void correctControl(const vectorField& U, const scalarField* rhoField) = 0;
protected:
    const word sourceName_;
    const fvMesh& mesh_;
};

}

#endif /* PROPELLER_PROPELLER_MODEL_PROPELLER_MODEL */
