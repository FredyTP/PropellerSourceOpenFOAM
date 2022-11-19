#ifndef PROPELLER_MODEL_BLADE_ELEMENT_MODEL_ELEMENT_MODEL_H
#define PROPELLER_MODEL_BLADE_ELEMENT_MODEL_ELEMENT_MODEL_H
#include "propellerModel.H"

namespace Foam
{

class bladeElementModel
:
    public propellerModel
{
public:
    //Defines the runtime name of this class (IOpurposes?)
    TypeName("bladeElementModel");

    // Constructors

        //- Construct from components
        bladeElementModel
        (
            const dictionary& dict
        );

                //- No copy construct
        bladeElementModel(const bladeElementModel&) = delete;

        //- No copy assignment
        void operator=(const bladeElementModel&) = delete;


    //- Destructor
    virtual ~bladeElementModel() = default;
    
};

}

#endif /* PROPELLER_MODEL_BLADE_ELEMENT_MODEL_ELEMENT_MODEL_H */

