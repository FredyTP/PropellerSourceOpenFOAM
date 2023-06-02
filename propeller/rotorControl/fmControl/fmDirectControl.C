#include "fmDirectControl.H"
#include "addToRunTimeSelectionTable.H"
#include "forceModel.H"

namespace Foam
{
 
defineTypeNameAndDebug(fmDirectControl,0);

addToRunTimeSelectionTable(fmControl,fmDirectControl,dictionary);

fmDirectControl::fmDirectControl(const dictionary& dict,const forceModel& fmModel)
 : fmControl(dict,fmModel)
{
    fixedJ_ = dict.readIfPresent("J",J_);
    if(!fixedJ_)
    {
        omega_ = rotorControl::readAngularVelocity(dict);
    }
}

void fmDirectControl::correctControl(const vectorField & U, const scalarField * rhoField)
{
    vector normal = forceModel_.grid()->geometry().direction();
    scalar diameter = 2 * forceModel_.grid()->geometry().radius();
    scalar speedRef = forceModel::getReferenceSpeed(U,normal);

    if(fixedJ_)
    {
        omega_ = speedRef/(J_/constant::mathematical::twoPi*diameter);
    }
    else
    {
        J_ = speedRef/(omega_/constant::mathematical::twoPi*diameter);
    }
}
scalar fmDirectControl::getAngularVelocity() const
{
    return omega_;
}
scalar fmDirectControl::getJ() const
{
    return J_;
}


}

