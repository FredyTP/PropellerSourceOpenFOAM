#ifndef PROPELLER_MODEL_FROUDE_MODEL_FROUDE_MODEL
#define PROPELLER_MODEL_FROUDE_MODEL_FROUDE_MODEL
#include "propellerModel.H"
#include "fmControl.H"
namespace Foam
{

class forceModel
:
    public propellerModel
{
public:
    //Defines the runtime name of this class (IOpurposes?)
    TypeName("forceModel");

    // Constructors

    //- Construct from components
    forceModel
    (
        const dictionary& dict
    );

            //- No copy construct
    forceModel(const forceModel&) = delete;

    //- No copy assignment
    void operator=(const forceModel&) = delete;


    //- Destructor
    virtual ~forceModel() = default;
    virtual void build(const rotorGeometry& rotorGeometry){}
    virtual void nextTimeStep(scalar dt);
    virtual propellerResult calculate(const vectorField& U, const scalarField* rhoField, volVectorField& force);
    virtual propellerResult calculate(const vectorField& U, const scalarField* rhoField) const;

    virtual rotorControl& control() {return control_.ref();};
    
    virtual const autoPtr<rotorGrid>& grid() {return rotorGrid_;}
private:
    dictionary gridDictionary_;
    autoPtr<rotorGrid> rotorGrid_;
    List<tensor> gridTensor_;
    autoPtr<fmControl> control_;

    autoPtr<regularInterpolation<scalar,scalar,1>> thrustCoeff_;
    autoPtr<regularInterpolation<scalar,scalar,1>> torqueCoeff_;

    static vector ForceDistribution(scalar Ax, scalar Atheta, scalar radius, scalar minRadius, scalar maxRadius);

    static scalar AxCoefficient(scalar thrust,scalar minRadius,scalar maxRadius);
    static scalar AthetaCoefficient(scalar torque,scalar minRadius,scalar maxRadius);


};

}

#endif /* PROPELLER_MODEL_FROUDE_MODEL_FROUDE_MODEL */
