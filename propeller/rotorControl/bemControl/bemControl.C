#include "bemControl.H"
#include "addToRunTimeSelectionTable.H"
#include "bladeElementModel.H"

namespace Foam
{
 
defineTypeNameAndDebug(bemControl,0);
defineRunTimeSelectionTable(bemControl, dictionary);

bemControl::bemControl(const dictionary& dict)
{
    
}
scalar bemControl::readAngularVelocity(const dictionary &dict)
{
    scalar rpm;

    if(dict.readIfPresent("rpm",rpm))
    {
        return rpmToRads(rpm);
    }
    else
    {
        return dict.get<scalar>("angularVelocity");
    }
}
vector bemControl::getBladeLocalVel(scalar azimuth, scalar radius) const
{
    return vector
    (
        -getAzimuthDot(azimuth)*radius, //x-axis "backwards"
        0,
        getFlappingDot(azimuth)*radius  //z-axis "up"
    );
}

autoPtr<bemControl> bemControl::New(const dictionary &dict, const bladeElementModel &bem)
{
    const word modelType(dict.get<word>("type")); 

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            typeName,
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<bemControl>(ctorPtr(dict, bem));
}

scalar bemControl::PitchFunction(scalar azimuth, scalar pitch0, scalar pitch1c, scalar pitch1s)
{
    return pitch0+pitch1c*cos(azimuth)+pitch1s*sin(azimuth);
}

scalar bemControl::AzimuthFunction(scalar azimuth0, scalar azimuth1c, scalar azimuth1s)
{
    return azimuth0+azimuth1c*cos(azimuth0)+azimuth1s*sin(azimuth0);
}
scalar bemControl::FlappingFunction(scalar azimuth, scalar flapping0, scalar flapping1c, scalar flapping1s)
{
    return flapping0+flapping1c*cos(azimuth)+flapping1s*sin(azimuth);
}
scalar bemControl::PitchDotFunction(scalar omega, scalar azimuth, scalar pitch0, scalar pitch1c, scalar pitch1s)
{
    return omega*(- pitch1c*sin(azimuth) + pitch1s*cos(azimuth) );
}
scalar bemControl::AzimuthDotFunction(scalar omega, scalar azimuth, scalar azimuth1c, scalar azimuth1s)
{
    return omega*(1 - azimuth1c*sin(azimuth) + azimuth1s*cos(azimuth) );
}
scalar bemControl::FlappingDotFunction(scalar omega, scalar azimuth, scalar flapping0, scalar flapping1c, scalar flapping1s)
{
    return omega*(- flapping1c*sin(azimuth) + flapping1s*cos(azimuth) );
}

}

