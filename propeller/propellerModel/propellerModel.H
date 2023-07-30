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
/**
 * @brief Base class for all propellerModels implementations
 * 
 */
class propellerModel
{

protected:
    const rotorFvMeshSel* rotorFvMeshSel_;
    autoPtr<rotorGrid> rotorGrid_;
    label nBlades_;
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

        /**
         * @brief Return a reference to the selected propellerModel
         * automations the construction of propeller models
         *   from run time selection tables
         * @param sourceName 
         * @param dict 
         * @param mesh 
         * @return autoPtr<propellerModel> 
         */
        static autoPtr<propellerModel> New
        (
            word sourceName,
            const dictionary& dict,
            const fvMesh& mesh
        );

    //Destructor

        virtual ~propellerModel() = default;


    void setRotorMesh(const rotorFvMeshSel* rotorFvMeshSel){ rotorFvMeshSel_ = rotorFvMeshSel;}

    virtual void build(const rotorGeometry& rotorGeometry)=0;
    virtual bool nextTimeStep(scalar dt) = 0;
    virtual propellerResult calculate(const vectorField& U, const scalarField* rhoField, volVectorField& force) = 0;
    virtual propellerResult calculate(const vectorField& U, const scalarField* rhoField) const =0;

    const autoPtr<rotorGrid>& grid() const {return rotorGrid_;}
    autoPtr<rotorGrid>& grid() {return rotorGrid_;}
    
    const fvMesh& mesh() const { return mesh_;}
    const word& sourceName() const {return sourceName_;}
    
    static tensor bladeTensor(const coordSystem::cylindrical &cylCS, const point &localPoint, scalar flapping, scalar sweep);

    /**
     * @brief Correct control values 
     * 
     * @param U 
     * @param rhoField 
     */
    virtual void correctControl(const vectorField& U, const scalarField* rhoField) = 0;
protected:
    const word sourceName_;
    const fvMesh& mesh_;
};

}

#endif /* PROPELLER_PROPELLER_MODEL_PROPELLER_MODEL */
