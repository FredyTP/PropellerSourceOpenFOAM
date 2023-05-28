#include "fmControl.H"
#include "addToRunTimeSelectionTable.H"
#include "forceModel.H"

namespace Foam
{
 
defineTypeNameAndDebug(fmControl,0);
defineRunTimeSelectionTable(fmControl, dictionary);

fmControl::fmControl(const dictionary& dict)
{
    
}
scalar fmControl::readAngularVelocity(const dictionary &dict)
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
autoPtr<fmControl> fmControl::New(const dictionary &dict, const forceModel &bem)
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

    return autoPtr<fmControl>(ctorPtr(dict, bem));
}

}

