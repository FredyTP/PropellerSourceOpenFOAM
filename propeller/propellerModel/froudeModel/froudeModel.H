#ifndef PROPELLER_MODEL_FROUDE_MODEL_FROUDE_MODEL
#define PROPELLER_MODEL_FROUDE_MODEL_FROUDE_MODEL
#include "propellerModel.H"

namespace Foam
{

class froudeModel
:
    public propellerModel
{
public:
    //Defines the runtime name of this class (IOpurposes?)
    TypeName("froudeModel");

    // Constructors

        //- Construct from components
        froudeModel
        (
            const dictionary& dict
        );

                //- No copy construct
        froudeModel(const froudeModel&) = delete;

        //- No copy assignment
        void operator=(const froudeModel&) = delete;


    //- Destructor
    virtual ~froudeModel() = default;

    
    virtual const autoPtr<rotorGrid>& grid() {return rotorGrid_;}
private:
    autoPtr<rotorGrid> rotorGrid_;
};

}

#endif /* PROPELLER_MODEL_FROUDE_MODEL_FROUDE_MODEL */
