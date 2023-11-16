#ifndef PROPELLER_PROPELLER_MODEL_BLADE_ELEMENT_MODEL_BLADE_ELEMENT_MODEL
#define PROPELLER_PROPELLER_MODEL_BLADE_ELEMENT_MODEL_BLADE_ELEMENT_MODEL

#include "propellerModel.H"
#include "airfoilModelList.H"
#include "bladeModelS.H"
#include "rotorGrid.H"
#include "polarCorrection.H"
#include "outputFields.H"

namespace Foam
{
struct bemDebugData
{
    scalar cl,cd,phi,aoa,radius,chord,twist,sweep;
    vector localVel;
};
class bemControl;

/**
 * @brief implementation of the blade element teory method. Valid for stationary or transient
 * 
 */
class bladeElementModel
:
    public propellerModel
{

public:
    /**PUBLIC ENUMERATIONS AND ENUM CLASS*/

    //- availabe scalar fiels
    enum outputScalarFields
    {
        aoaField=0,
        clField,
        cdField,
    };

    typedef util::outputFields<scalar,outputScalarFields,3> outputScalarType;
    //- available vector fields
    enum outputVectorFields
    {
        bemForceField,
        bemMomentField,
        bladeVelocityField,
        inflowVelocityField,
        relativeVelocityField
    };
    
    typedef util::outputFields<vector,outputVectorFields,5> outputVectorType;


    enum class controlVar
    {
        omega=0,
        collectivePitch,
        ciclicPitchCos,
        ciclicPitchSin
    };
    enum class outputVar
    {
        forceX=0,
        forceY,
        forceZ,
        torqueX,
        torqueY,
        torqueZ,
        power
    };
    const static Enum<controlVar> controlVarNames_;
    const static Enum<outputVar> outputVarNames_;
    
    static const label N_control = 4;
    static const label N_output = 7;
protected:

    airfoilModelList airfoils_;
    bladeModelS bladeModel_;
    List<tensor> gridTensor_; //(x: LE->TE) (y: radial) (z: "up")
    scalar tipFactor_;

    scalar rhoRef_;
    scalar speedRef_;
    scalar nuRef_;
    scalar soundSpeedRef_;

    //- Reference frame rotation for time-accurate blades
    scalar psi0_ = 0;

    dictionary gridDictionary;
    autoPtr<bemControl> control_;
    effects3D::polarCorrection polarCorrection_;

    mutable outputScalarType scalarFields_;
    mutable outputVectorType vectorFields_;

public:
    //Defines the runtime name of this class (IOpurposes?)
    TypeName("bladeElementModel");

    
    // Constructors

        //- Construct from components
        bladeElementModel
        (
            word sourceName,
            const dictionary& dict,
            const fvMesh& mesh
        );

            //- No copy construct
        bladeElementModel(const bladeElementModel&) = delete;

        //- No copy assignment
        void operator=(const bladeElementModel&) = delete;


    //- Destructor
    virtual ~bladeElementModel() = default;

    //-Getter


    virtual void build(const rotorGeometry& rotorGeometry);

    /**
     * @brief Calculates and updates theta state, also add source term
     * 
     * @param U 
     * @param rhoField 
     * @param angularVelocity 
     * @param force 
     * @param theta 
     * @return propellerResult 
     */
    virtual propellerResult calculate(const vectorField& U, const scalarField* rhoField, volVectorField& force);

    /**
     * @brief Calculates only the forces on the current state
     * 
     * @param U 
     * @param rhoField 
     * @param angularVelocity 
     * @return propellerResult 
     */
    virtual propellerResult calculate(const vectorField& U, const scalarField* rhoField) const;

    

    vector calculatePoint(const vector& U, scalar rho, scalar angularVelocity, const gridCell& cell, const tensor& bladeTensor, bemDebugData& data, bool saveFields = false) const;

    virtual bool nextTimeStep(scalar dt);
    
    static inline scalar AngleOfIncidenceSTAR(const vector& relativeLocalVel);
    static inline scalar AngleOfIncidence(const vector& relativeLocalVel);

    inline scalar rhoRef() const { return rhoRef_;}

    bool isTimeAcuratte() const;

    void updateTensors();

    tensor cellBladeTensor(const gridCell& cell) const;
    
    void correctControl(const vectorField& U, const scalarField* rhoField);
};

}

#endif /* PROPELLER_PROPELLER_MODEL_BLADE_ELEMENT_MODEL_BLADE_ELEMENT_MODEL */

