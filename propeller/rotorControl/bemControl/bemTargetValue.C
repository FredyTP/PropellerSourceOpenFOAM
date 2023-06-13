#include "bemTargetValue.H"
#include "unitConversion.H"
#include "functionSolver.H"
namespace Foam 
{
defineTypeNameAndDebug(bemTargetValue,0);
addToRunTimeSelectionTable(bemControl,bemTargetValue, dictionary);


bemTargetValue::bemTargetValue(const dictionary &dict, const bladeElementModel& bem)
:
    bemControl(dict),
    bem_(bem),
    calcFrequency_(-1),
    target_(),
    control_(),
    nIter_(50),
    tol_(1e-8),
    relax_(1.0),
    dx_(1,0)
{
    read(dict);
}

void bemTargetValue::read(const dictionary &dict)
{
    //Read target variables and get the list of the targets used
    usedTarget_ = target_.readIfPresent(dict,bladeElementModel::outputVarNames_);

    //Read control variables initial value for thos used and state for thos unused
    control_.readOrDefault(dict,bladeElementModel::controlVarNames_,0);
    control_[controlVar::omega] = rotorControl::readAngularVelocity(dict);

    //Get used control variables
    List<word> controlNames = dict.get<List<word>>("controlVariables");
    
    usedControl_.resize(controlNames.size());
    forAll(controlNames, i)
    {
        usedControl_[i] = bladeElementModel::controlVarNames_.get(controlNames[i]);
    }

    //Read simulation properties
    dict.readEntry("calcFrequency", calcFrequency_);
    dict.readIfPresent("nIter", nIter_);
    dict.readIfPresent("tol", tol_);
    dict.readIfPresent("relax", relax_);

    if (!dict.readIfPresent("dx", dx_))
    {
        dx_ = scalarField(usedControl_.size(), 0.01);
    }
}

void bemTargetValue::correctControl(const vectorField &U, const scalarField *rhoField)
{
    if (bem_.mesh().time().timeIndex() % calcFrequency_ == 0)
    {
        word calcType = "forces";

        Info<< type() << ":" << nl
            << "    solving for target trim " << calcType << nl;

        const scalar rhoRef = bem_.rhoRef();
        scalarField x0;
        control_.get(x0,usedControl_);
  
        util::functionSolver::NewtonRapson
        (
            usedControl_.size(),
            solverFunction(U,rhoField),
            x0,
            dx_,
            relax_,
            nIter_,
            tol_,
            true
        );

        Info<<control_<<endl;
    }
}


vector bemTargetValue::calcForces(const vectorField &U, const scalarField *rhoField) const
{
    vector forces;
    auto result = bem_.calculate(U,rhoField);
    forces[0] = result.thrust();
    forces[1] = result.torque.x();
    forces[2] = result.torque.y();
    return forces;
}


scalar bemTargetValue::getAzimuth(scalar azimuth0) const
{
    return azimuth0;
}

scalar bemTargetValue::getPitch(scalar azimuth) const
{
    return PitchFunction
    (
        azimuth,
        control_[controlVar::collectivePitch],
        control_[controlVar::ciclicPitchCos],
        control_[controlVar::ciclicPitchSin]
    );
}

scalar bemTargetValue::getFlapping(scalar azimuth) const
{
    return 0;
}

scalar bemTargetValue::getAzimuthDot(scalar azimuth) const
{
    return this->getAngularVelocity();
}
scalar bemTargetValue::getPitchDot(scalar azimuth) const
{
    return PitchDotFunction
    (
        control_[controlVar::omega],
        azimuth,
        control_[controlVar::collectivePitch],
        control_[controlVar::ciclicPitchCos],
        control_[controlVar::ciclicPitchSin]
    );
}
scalar bemTargetValue::getFlappingDot(scalar azimuth) const
{
    return 0;
}

scalar bemTargetValue::getAngularVelocity() const
{
    return control_[controlVar::omega];
}

std::function<scalarField(scalarField)> bemTargetValue::solverFunction(const vectorField &U, const scalarField *rhoField)
{
    return [this,U,rhoField](const scalarField& control)
    {
        
        control_.set(control,usedControl_);
        
        auto result = bem_.calculate(U,rhoField);
        outputVarType output = outputFromResult(result);

        scalarField values;
        output.get(values,usedTarget_);
        scalarField targetValues;
        target_.get(targetValues,usedTarget_);

        values = values - targetValues;
        return values;
    };
}
typename bemTargetValue::outputVarType bemTargetValue::outputFromResult(const propellerResult &result) const
{
    outputVarType ret;

    ret[outputVar::forceX] = result.force.x();
    ret[outputVar::forceY] = result.force.y();
    ret[outputVar::forceZ] = result.force.z();
    ret[outputVar::torqueX] = result.torque.x();
    ret[outputVar::torqueY] = result.torque.y();
    ret[outputVar::torqueZ] = result.torque.z();
    ret[outputVar::power] = result.power;

    return ret;
}
}
