#ifndef PROPELLER_MODEL_PROPELLER_MODEL
#define PROPELLER_MODEL_PROPELLER_MODEL

#endif /* PROPELLER_MODEL_PROPELLER_MODEL */

#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "rotorMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class propellerModel Declaration
\*---------------------------------------------------------------------------*/




class propellerModel
{

protected:
    const rotorMesh* rotorMesh_;
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

    virtual void calculate(const vectorField& U,volVectorField& force){}
};

}
