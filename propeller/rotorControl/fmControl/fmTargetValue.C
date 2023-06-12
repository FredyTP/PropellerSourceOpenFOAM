#include "fmTargetValue.H"
#include "addToRunTimeSelectionTable.H"
#include "forceModel.H"
#include "functionSolver.H"
namespace Foam
{
 
defineTypeNameAndDebug(fmTargetValue,0);

addToRunTimeSelectionTable(fmControl,fmTargetValue,dictionary);

fmTargetValue::fmTargetValue(const dictionary& dict,const forceModel& fmModel)
 : fmControl(dict,fmModel),
    calcFrequency_(-1),
    target_(),
    control_(),
    usedControl_(1,controlVar::omega),
    nIter_(50),
    tol_(1e-8),
    relax_(1.0)
{
    read(dict);
}
void fmTargetValue::read(const dictionary &dict)
{
    //Read target variables and get the list of the targets used
    usedTarget_ = target_.readIfPresent(dict,forceModel::outputVarNames_);

    updateOmega_ = dict.readIfPresent("J",initalJ_);
    if(!updateOmega_)
    {
        control_[controlVar::omega] = rotorControl::readAngularVelocity(dict);
    }

    

    //Read simulation properties
    dict.readEntry("calcFrequency", calcFrequency_);
    dict.readIfPresent("nIter", nIter_);
    dict.readIfPresent("tol", tol_);
    dict.readIfPresent("relax", relax_);
}


void fmTargetValue::correctControl(const vectorField & U, const scalarField * rhoField)
{
    if(updateOmega_)
    {
        updateOmega_=false;
        control_[controlVar::omega] = fmControl::getOmegaFromJ(initalJ_,this->forceModel_.getReferenceSpeed(U));
    }
    if (this->forceModel_.mesh().time().timeIndex() % calcFrequency_ == 0)
    {
        refSpeed_ = this->forceModel_.getReferenceSpeed(U);
        word calcType = "forces";

        Info<< type() << ":" << nl
            << "    solving for target trim " << calcType << nl;

        scalarField x0;
        control_.get(x0,usedControl_);
        util::functionSolver::NewtonRapson
        (
            usedControl_.size(),
            solverFunction(U,rhoField),
            x0,
            scalarField(usedControl_.size(),4),
            relax_,
            nIter_,
            tol_,
            true
        );

    }
}


std::function<scalarField(scalarField)> fmTargetValue::solverFunction(const vectorField &U, const scalarField *rhoField)
{
    return [this,U,rhoField](const scalarField& control)
    {
        
        control_.set(control,usedControl_);
        
        auto result = this->forceModel_.calculate(U,rhoField);
        outputVarType output = outputFromResult(result);

        scalarField values;
        output.get(values,usedTarget_);
        scalarField targetValues;
        target_.get(targetValues,usedTarget_);

        values = values - targetValues;
        return values;
    };
}
typename fmTargetValue::outputVarType fmTargetValue::outputFromResult(const propellerResult &result) const
{
    outputVarType ret;
    ret[outputVar::forceZ] = result.force.z();
    ret[outputVar::torqueZ] = result.torque.z();
    ret[outputVar::power] = result.power;

    return ret;
}

scalar fmTargetValue::getAngularVelocity() const
{
    return control_[controlVar::omega];
}
scalar fmTargetValue::getJ() const
{
    return fmControl::getJFromOmega(control_[controlVar::omega],refSpeed_);
}


}

