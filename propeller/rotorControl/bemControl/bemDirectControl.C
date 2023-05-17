#include "bemDirectControl.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
 
defineTypeNameAndDebug(bemDirectControl,0);
addToRunTimeSelectionTable(rotorControl,bemDirectControl, dictionary);


bemDirectControl::bemDirectControl(const dictionary& dict)
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
}

scalar bemDirectControl::getAzimuth(scalar azimuth0) const
{
    return azimuth0-sinAzimuth*sin(azimuth0)-cosAzimuth*cos(azimuth0);
}

scalar bemDirectControl::getPitch(scalar azimuth) const
{
    return collectivePitch-sinPitch*sin(azimuth)-cosPitch*cos(azimuth);
}

scalar bemDirectControl::getFlapping(scalar azimuth) const
{
    return flapping-sinFlapping*sin(azimuth)-cosFlapping*cos(azimuth);
}

scalar bemDirectControl::getAzimuthDot(scalar azimuth, scalar angularVelocity)
{
    return angularVelocity*(1 - sinAzimuth*cos(azimuth) + cosAzimuth*sin(azimuth));
}
scalar bemDirectControl::getPitchDot(scalar azimuth, scalar angularVelocity)
{
    return angularVelocity*(-sinPitch*cos(azimuth) + cosPitch*sin(azimuth));
}
scalar bemDirectControl::getFlappingDot(scalar azimuth, scalar angularVelocity)
{
    return angularVelocity*(-sinFlapping*cos(azimuth) + cosFlapping*sin(azimuth));
}

vector bemDirectControl::getBladeLocalVel(scalar azimuth, scalar angularVelocity, scalar radius)
{
    return vector
    (
        -getAzimuthDot(azimuth,angularVelocity)*radius, //x-axis "backwards"
        0,
        getFlappingDot(azimuth,angularVelocity)*radius  //z-axis "up"
    );
}
}