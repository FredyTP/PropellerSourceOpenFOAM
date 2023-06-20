#ifndef PROPELLER_PROPELLER_MODEL_FORCE_MODEL_FORCE_MODEL
#define PROPELLER_PROPELLER_MODEL_FORCE_MODEL_FORCE_MODEL
#include "propellerModel.H"
#include "RegularInterpolation.H"
#include "vector.H"
#include "tensor.H"

namespace Foam
{
class fmControl;

class forceModel
:
    public propellerModel
{
public:
    //Defines the runtime name of this class (IOpurposes?)
    TypeName("forceModel");

    enum class controlVar
    {
        omega=0
    };
    enum class outputVar
    {
        forceZ=0,
        torqueZ,
        power
    };

    static const label N_control = 2;
    static const label N_output = 3;
    
    const static Enum<controlVar> controlVarNames_;
    const static Enum<outputVar> outputVarNames_;
    // Constructors

    //- Construct from components
    forceModel
    (
        word sourceName,
        const dictionary& dict,
        const fvMesh& mesh
    );

    //- No copy construct
    forceModel(const forceModel&) = delete;

    //- No copy assignment
    void operator=(const forceModel&) = delete;


    //- Destructor
    virtual ~forceModel() = default;

    void read(const dictionary& dict);

    virtual void build(const rotorGeometry& rotorGeometry);
    virtual bool nextTimeStep(scalar dt);
    virtual propellerResult calculate(const vectorField& U, const scalarField* rhoField, volVectorField& force);
    virtual propellerResult calculate(const vectorField& U, const scalarField* rhoField) const;


    //virtual rotorControl<forceModel>& control() {return control_.ref();};
    
    virtual const autoPtr<rotorGrid>& grid() const {return rotorGrid_;}
private:
    
    scalar psi0_ = 0;
    scalar rhoRef_;
    dictionary gridDictionary_;
    autoPtr<rotorGrid> rotorGrid_;
    List<tensor> gridTensor_;
    autoPtr<fmControl> control_;

    autoPtr<RegularInterpolation<scalar,scalar,1>> thrustCoeff_;
    autoPtr<RegularInterpolation<scalar,scalar,1>> torqueCoeff_;

    static vector ForceDistribution(scalar Ax, scalar Atheta, scalar radius, scalar minRadius, scalar maxRadius);
    static vector ForceIntergralFunction(scalar Ax, scalar Atheta, scalar radius, scalar minRadius, scalar maxRadius);
    static vector ForceIntergral(scalar Ax, scalar Atheta, scalar radius0, scalar radius1, scalar minRadius, scalar maxRadius);
    static scalar AxCoefficient(scalar thrust,scalar minRadius,scalar maxRadius);
    static scalar AthetaCoefficient(scalar torque,scalar minRadius,scalar maxRadius);

    void updateTensors();
    tensor cellBladeTensor(const gridCell& cell) const;


    
public:
    scalar rhoRef() const {return rhoRef_;}
    static scalar getReferenceSpeed(const vectorField& U,const vector& normal);
    scalar getReferenceSpeed(const vectorField& U) const;

};

}

#endif /* PROPELLER_PROPELLER_MODEL_FORCE_MODEL_FORCE_MODEL */
