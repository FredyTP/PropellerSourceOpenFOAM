#include "fmControl.H"
#include "addToRunTimeSelectionTable.H"
#include "forceModel.H"

namespace Foam
{
 
defineTypeNameAndDebug(fmControl,0);
defineRunTimeSelectionTable(fmControl, dictionary);

fmControl::fmControl(const dictionary& dict,const forceModel& fmModel)
 : forceModel_(fmModel)
{
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

scalar fmControl::getJFromOmega(scalar omega, scalar speed) const
{
    scalar diameter = 2 * forceModel_.grid()->geometry().radius();

    return speed/(omega/constant::mathematical::twoPi*diameter);
}
scalar fmControl::getOmegaFromJ(scalar J, scalar speed) const
{
    scalar diameter = 2 * forceModel_.grid()->geometry().radius();

    return speed/(J/constant::mathematical::twoPi*diameter);
}


}
