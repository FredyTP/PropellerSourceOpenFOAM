#ifndef PROPELLER_PROPELLER_MODEL_BLADE_ELEMENT_MODEL_BLADE_ELEMENT_MODEL
#define PROPELLER_PROPELLER_MODEL_BLADE_ELEMENT_MODEL_BLADE_ELEMENT_MODEL
#include "propellerModel.H"
#include "airfoilModelList.H"
#include "bladeModelDev.H"
#include "rotorDiscrete.H"

namespace Foam
{

class bladeElementModel
:
    public propellerModel
{
protected:

    airfoilModelList airfoils_;
    Foam::devel::bladeModel bladeModel_;
    rotorDiscrete rotorDiscrete_;
    
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

    //-Getter
    virtual scalar radius() const;

    virtual void build(const rotorGeometry& rotorGeometry);

    virtual void calculate(volVectorField& force);
};

}

#endif /* PROPELLER_PROPELLER_MODEL_BLADE_ELEMENT_MODEL_BLADE_ELEMENT_MODEL */
