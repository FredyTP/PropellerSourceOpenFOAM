#include "bemTargetValue.H"
#include "bladeElementModel.H"
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
    dTheta_(degToRad(0.1))
{
    read(dict);
}

void bemTargetValue::read(const dictionary &dict)
{
    dict.readEntry("thrust", target_[outputVar::forceZ]);
    dict.readEntry("pitch", target_[outputVar::torqueX]);
    dict.readEntry("roll", target_[outputVar::torqueY]);

    control_[controlVar::collectivePitch] = degToRad(dict.get<scalar>("theta0Ini"));
    control_[controlVar::ciclicPitchCos] = degToRad(dict.get<scalar>("theta1cIni"));
    control_[controlVar::ciclicPitchSin] = degToRad(dict.get<scalar>("theta1sIni"));

    dict.readEntry("calcFrequency", calcFrequency_);

    dict.readIfPresent("nIter", nIter_);
    dict.readIfPresent("tol", tol_);
    dict.readIfPresent("relax", relax_);

    if (dict.readIfPresent("dTheta", dTheta_))
    {
        dTheta_ = degToRad(dTheta_);
    }

    control_[controlVar::omega] = rotorControl::readAngularVelocity(dict);

    usedControl_.resize(3);
    usedTarget_.resize(3);

    usedControl_[0]=controlVar::collectivePitch;
    usedControl_[1]=controlVar::ciclicPitchCos;
    usedControl_[2]=controlVar::ciclicPitchSin;

    usedTarget_[0]=outputVar::forceZ;
    usedTarget_[1]=outputVar::torqueX;
    usedTarget_[2]=outputVar::torqueY;

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
        util::functionSolver<3>::NewtonRapson(solverFunction(U,rhoField),x0,scalarField(3,dTheta_),relax_,nIter_,tol_,true);
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
        return output;
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

    return ret;
}
}
