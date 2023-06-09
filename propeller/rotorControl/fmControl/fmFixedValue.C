#include "fmFixedValue.H"
#include "addToRunTimeSelectionTable.H"
#include "forceModel.H"

namespace Foam
{
 
defineTypeNameAndDebug(fmFixedValue,0);

addToRunTimeSelectionTable(fmControl,fmFixedValue,dictionary);

fmFixedValue::fmFixedValue(const dictionary& dict,const forceModel& fmModel)
 : fmControl(dict,fmModel)
{
    fixedJ_ = dict.readIfPresent("J",J_);
    if(!fixedJ_)
    {
        omega_ = rotorControl::readAngularVelocity(dict);
    }
}

void fmFixedValue::correctControl(const vectorField & U, const scalarField * rhoField)
{
    vector normal = forceModel_.grid()->geometry().direction();
    scalar speedRef = forceModel::getReferenceSpeed(U,normal);

    if(fixedJ_)
    {
        omega_ = fmControl::getOmegaFromJ(J_,speedRef);
    }
    else
    {
        J_ = fmControl::getJFromOmega(omega_,speedRef);
    }
}
scalar fmFixedValue::getAngularVelocity() const
{
    return omega_;
}
scalar fmFixedValue::getJ() const
{
    return J_;
}


}

