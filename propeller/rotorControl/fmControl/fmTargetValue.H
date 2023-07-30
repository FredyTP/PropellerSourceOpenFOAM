#ifndef PROPELLER_ROTOR_CONTROL_FM_CONTROL_FM_TARGET_VALUE
#define PROPELLER_ROTOR_CONTROL_FM_CONTROL_FM_TARGET_VALUE

#include "fmControl.H"
#include "unitConversion.H"
#include "vectorField.H"
#include "scalarField.H"

namespace Foam
{

class forceModel;

/**
 * @brief Class for force model rotor control of target output values via Newton-Raphson algorithm
 * 
 */
class fmTargetValue : public fmControl
{
public:

    TypeName("targetValue");

    fmTargetValue(const dictionary& dict, const forceModel& fmModel);

    virtual void correctControl(const vectorField& U, const scalarField* rhoField);
    
    void read(const dictionary &dict);

    virtual scalar getAngularVelocity() const;

    virtual scalar getJ() const;

    std::function<scalarField(scalarField)> solverFunction(const vectorField&u, const scalarField *rhoField);

    outputVarType outputFromResult(const propellerResult& result) const;
 
private:
    scalar refSpeed_;
    bool updateOmega_ = false;
    bool updateJ_ = false;
    scalar initalJ_;
    //- Number of iterations between calls to 'correct'
    label calcFrequency_;

    //- Target coefficient vector (thrust force, roll moment, pitch moment)
    outputVarType target_;

    //- Pitch angles (collective, roll, pitch) [rad]
    controlVarType control_;

    //- Used variables
    List<outputVar> usedTarget_;
    List<controlVar> usedControl_;

    //- Maximum number of iterations in trim routine
    label nIter_;

    //- Relaxation factor
    scalar relax_;

    //- Convergence tolerance
    scalar tol_;

    //- Dx used
    scalarList dx_;

    

};

}




#endif /* PROPELLER_ROTOR_CONTROL_FM_CONTROL_FM_TARGET_VALUE */
