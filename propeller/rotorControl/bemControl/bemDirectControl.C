#include "bemDirectControl.H"
#include "addToRunTimeSelectionTable.H"
#include "bladeElementModel.H"
namespace Foam
{
 
defineTypeNameAndDebug(bemDirectControl,0);
addToRunTimeSelectionTable(bemControl,bemDirectControl, dictionary);


bemDirectControl::bemDirectControl(const dictionary& dict, const bladeElementModel& bem)
 : bemControl(dict)
{
    collectivePitch = dict.getOrDefault<scalar>("collectivePitch",0);
    sinPitch = dict.getOrDefault<scalar>("sinPitch",0);
    cosPitch = dict.getOrDefault<scalar>("cosPitch",0);
    sinAzimuth = dict.getOrDefault<scalar>("sinAzimuth",0);
    cosAzimuth = dict.getOrDefault<scalar>("cosAzimuth",0);
    flapping = dict.getOrDefault<scalar>("flapping",0);
    sinFlapping = dict.getOrDefault<scalar>("sinFlapping",0);
    cosFlapping = dict.getOrDefault<scalar>("cosFlapping",0);

    angularVelocity_ = rotorControl::readAngularVelocity(dict);
}

void bemDirectControl::correctControl(const vectorField &U, const scalarField *rhoField)
{

}

scalar bemDirectControl::getAzimuth(scalar azimuth0) const
{
    return AzimuthFunction(azimuth0,cosAzimuth,sinAzimuth);
}

scalar bemDirectControl::getPitch(scalar azimuth) const
{
    return PitchFunction(azimuth,collectivePitch,cosPitch,sinPitch);
}

scalar bemDirectControl::getFlapping(scalar azimuth) const
{
    return FlappingFunction(azimuth,flapping,cosFlapping,sinFlapping);
}

scalar bemDirectControl::getAzimuthDot(scalar azimuth) const
{
    return AzimuthDotFunction(angularVelocity_,azimuth,cosAzimuth,sinAzimuth);
}
scalar bemDirectControl::getPitchDot(scalar azimuth) const
{
    return PitchDotFunction(angularVelocity_,azimuth,collectivePitch,cosPitch,sinPitch);
}
scalar bemDirectControl::getFlappingDot(scalar azimuth) const
{
    return FlappingDotFunction(angularVelocity_,azimuth,flapping,cosFlapping,sinFlapping);
}

scalar bemDirectControl::getAngularVelocity() const
{
    return angularVelocity_;
}
}