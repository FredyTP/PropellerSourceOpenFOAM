#ifndef PROPELLER_MODEL_FROUDE_MODEL_FROUDE_MODEL
#define PROPELLER_MODEL_FROUDE_MODEL_FROUDE_MODEL
#include "propellerModel.H"
#include "fmControl.H"
#include "vector.H"
#include "tensor.H"

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
    virtual void build(const rotorGeometry& rotorGeometry);
    virtual void nextTimeStep(scalar dt);
    virtual propellerResult calculate(const vectorField& U, const scalarField* rhoField, volVectorField& force);
    virtual propellerResult calculate(const vectorField& U, const scalarField* rhoField) const;


    virtual rotorControl& control() {return control_.ref();};
    
    virtual const autoPtr<rotorGrid>& grid() const {return rotorGrid_;}
private:
    
    scalar psi0_ = 0;
    scalar rhoRef_;
    dictionary gridDictionary_;
    autoPtr<rotorGrid> rotorGrid_;
    List<tensor> gridTensor_;
    autoPtr<fmControl> control_;
    bladeModelS bladeModel_;

    autoPtr<regularInterpolation<scalar,scalar,1>> thrustCoeff_;
    autoPtr<regularInterpolation<scalar,scalar,1>> torqueCoeff_;

    static vector ForceDistribution(scalar Ax, scalar Atheta, scalar radius, scalar minRadius, scalar maxRadius);

    static scalar AxCoefficient(scalar thrust,scalar minRadius,scalar maxRadius);
    static scalar AthetaCoefficient(scalar torque,scalar minRadius,scalar maxRadius);

    void updateTensors();
    tensor cellBladeTensor(const gridCell& cell) const;

};

}

#endif /* PROPELLER_MODEL_FROUDE_MODEL_FROUDE_MODEL */
