#include "fmControl.H"
#include "addToRunTimeSelectionTable.H"
#include "forceModel.H"

namespace Foam
{
 
defineTypeNameAndDebug(fmControl,0);
defineRunTimeSelectionTable(fmControl, dictionary);

addToRunTimeSelectionTable(fmControl,fmControl,dictionary);

fmControl::fmControl(const dictionary& dict,const forceModel& fmModel)
 : forceModel_(fmModel)
{
    omega_ = readAngularVelocity(dict);
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
void fmControl::correctControl(const vectorField & U, const scalarField * rhoField)
{
    vector normal = forceModel_.grid()->geometry().direction();
    scalar diameter = 2* forceModel_.grid()->geometry().radius();
 
    vector velAvg = average(U);
    Info<<"vel Avg "<<velAvg<<endl;
    scalar speedRef = -normal.inner(velAvg);
    Info<<"speed: "<<speedRef<<endl;
    J=speedRef/(omega_/constant::mathematical::twoPi*diameter);
}
scalar fmControl::getJ()
{
    return J;
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

