#include "bemFixedControl.H"
#include "addToRunTimeSelectionTable.H"
#include "bladeElementModel.H"
namespace Foam
{
 
defineTypeNameAndDebug(bemFixedControl,0);
addToRunTimeSelectionTable(bemControl,bemFixedControl, dictionary);


bemFixedControl::bemFixedControl(const dictionary& dict, const bladeElementModel& bem)
 : bemControl(dict)
{
    control_.readOrDefault(dict,bladeElementModel::controlVarNames_,0);
    sinAzimuth = dict.getOrDefault<scalar>("sinAzimuth",0);
    cosAzimuth = dict.getOrDefault<scalar>("cosAzimuth",0);
    flapping = dict.getOrDefault<scalar>("flapping",0);
    sinFlapping = dict.getOrDefault<scalar>("sinFlapping",0);
    cosFlapping = dict.getOrDefault<scalar>("cosFlapping",0);

    control_[controlVar::omega] = rotorControl::readAngularVelocity(dict);
}

void bemFixedControl::correctControl(const vectorField &U, const scalarField *rhoField)
{

}

scalar bemFixedControl::getAzimuth(scalar azimuth0) const
{
    return AzimuthFunction(azimuth0,cosAzimuth,sinAzimuth);
}

scalar bemFixedControl::getPitch(scalar azimuth) const
{
    return PitchFunction
    (
        azimuth,
        control_[controlVar::collectivePitch],
        control_[controlVar::ciclicPitchCos],
        control_[controlVar::ciclicPitchSin]
    );
}

scalar bemFixedControl::getFlapping(scalar azimuth) const
{
    return FlappingFunction(azimuth,flapping,cosFlapping,sinFlapping);
}

scalar bemFixedControl::getAzimuthDot(scalar azimuth) const
{
    return AzimuthDotFunction(control_[controlVar::omega],azimuth,cosAzimuth,sinAzimuth);
}
scalar bemFixedControl::getPitchDot(scalar azimuth) const
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
scalar bemFixedControl::getFlappingDot(scalar azimuth) const
{
    return FlappingDotFunction(control_[controlVar::omega],azimuth,flapping,cosFlapping,sinFlapping);
}

scalar bemFixedControl::getAngularVelocity() const
{
    return control_[controlVar::omega];
}
}