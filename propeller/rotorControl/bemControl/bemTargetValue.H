#ifndef PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_TARGET_VALUE
#define PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_TARGET_VALUE

#include "bemControl.H"
#include "EnumField.H"
#include "propellerResult.H"
#include "Enum.H"

namespace Foam
{

class bemTargetValue : public bemControl
{
public:
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

    typedef util::EnumField<outputVar,7> outputVarType;
    typedef util::EnumField<controlVar,4> controlVarType;

    TypeName("bemTargetValue");

    bemTargetValue(const dictionary& dict, const bladeElementModel& bem);

    void read(const dictionary& dict);
    
    void correctControl(const vectorField& U, const scalarField* rhoField);

    scalar getAzimuth(scalar azimuth0) const;
    scalar getPitch(scalar azimuth) const;
    scalar getFlapping(scalar azimuth) const;

    scalar getAzimuthDot(scalar azimuth) const;
    scalar getPitchDot(scalar azimuth) const;
    scalar getFlappingDot(scalar azimuth) const;

    vector calcForces(const vectorField &U, const scalarField *rhoField) const;
    virtual scalar getAngularVelocity() const;

    std::function<scalarField(scalarField)> solverFunction(const vectorField&u, const scalarField *rhoField);

    outputVarType outputFromResult(const propellerResult& result) const;
private:

    const bladeElementModel& bem_;

    //- Number of iterations between calls to 'correct'
    label calcFrequency_;

    //- Target coefficient vector (thrust force, roll moment, pitch moment)
    outputVarType target_;

    //- Pitch angles (collective, roll, pitch) [rad]
    controlVarType control_;

    List<outputVar> usedTarget_;
    List<controlVar> usedControl_;

    //- Maximum number of iterations in trim routine
    label nIter_;

    scalar relax_;

    scalar dTheta_;
    //- Convergence tolerance
    scalar tol_;
};




}


#endif /* PROPELLER_ROTOR_CONTROL_BEM_CONTROL_BEM_TARGET_VALUE */
