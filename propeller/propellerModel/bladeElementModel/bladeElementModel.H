#ifndef PROPELLER_PROPELLER_MODEL_BLADE_ELEMENT_MODEL_BLADE_ELEMENT_MODEL
#define PROPELLER_PROPELLER_MODEL_BLADE_ELEMENT_MODEL_BLADE_ELEMENT_MODEL
#include "propellerModel.H"
#include "airfoilModelList.H"
#include "bladeModelS.H"
#include "rotorDiscrete.H"
#include "velocitySampler.H"
namespace Foam
{
struct bemDebugData
{
    scalar cl,cd,phi,aoa,radius,chord,twist,sweep;
    vector localVel;
};
class bladeElementModel
:
    public propellerModel
{
protected:

    airfoilModelList airfoils_;
    bladeModelS bladeModel_;
    rotorDiscrete rotorDiscrete_;

    scalar nBlades_; 
    scalar tipFactor_;
    word integrationOrder;
    bool correctCenters;
    label refinementLevel;
    word discreteMethod;
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


    virtual void build(const rotorGeometry& rotorGeometry);

    virtual propellerResult calculate(const vectorField& U,scalar angularVelocity, volVectorField& force, scalar theta);
    
    virtual const rotorDiscrete& rDiscrete() {return rotorDiscrete_;}

    vector calculatePoint(const vector& U,scalar angularVelocity, const gridCell& cell, bemDebugData& data);

    
    static inline scalar AngleOfIncidenceSTAR(const vector& relativeLocalVel);
    static inline scalar AngleOfIncidence(const vector& relativeLocalVel);

    

};

}

#endif /* PROPELLER_PROPELLER_MODEL_BLADE_ELEMENT_MODEL_BLADE_ELEMENT_MODEL */

