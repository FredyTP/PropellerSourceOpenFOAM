#ifndef PROPELLER_PROPELLER_MODEL_BLADE_ELEMENT_MODEL_BLADE_ELEMENT_MODEL
#define PROPELLER_PROPELLER_MODEL_BLADE_ELEMENT_MODEL_BLADE_ELEMENT_MODEL
#include "propellerModel.H"
#include "airfoilModelList.H"
#include "bladeModelS.H"
#include "rotorGrid.H"
#include "polarCorrection.H"


namespace Foam
{
struct bemDebugData
{
    scalar cl,cd,phi,aoa,radius,chord,twist,sweep;
    vector localVel;
};
class bemControl;

class bladeElementModel
:
    public propellerModel
{
protected:

    airfoilModelList airfoils_;
    bladeModelS bladeModel_;
    autoPtr<rotorGrid> rotorGrid_;
    List<tensor> gridTensor_; //(x: LE->TE) (y: radial) (z: "up")
    scalar nBlades_; 
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
public:
    //Defines the runtime name of this class (IOpurposes?)
    TypeName("bladeElementModel");

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

    virtual const autoPtr<rotorGrid>& grid() const {return rotorGrid_;}

    vector calculatePoint(const vector& U, scalar rho, scalar angularVelocity, const gridCell& cell, const tensor& bladeTensor, bemDebugData& data) const;

    virtual void nextTimeStep(scalar dt);
    
    static inline scalar AngleOfIncidenceSTAR(const vector& relativeLocalVel);
    static inline scalar AngleOfIncidence(const vector& relativeLocalVel);

    inline scalar rhoRef() const { return rhoRef_;}

    bool isTimeAcuratte() const;

    void updateTensors();

    tensor cellBladeTensor(const gridCell& cell) const;
    
};

}

#endif /* PROPELLER_PROPELLER_MODEL_BLADE_ELEMENT_MODEL_BLADE_ELEMENT_MODEL */

