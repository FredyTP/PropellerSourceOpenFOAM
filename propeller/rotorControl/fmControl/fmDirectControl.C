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
scalar fmDirectControl::getAngularVelocity() const
{
    return omega_;
}
scalar fmDirectControl::getJ() const
{
    return J_;
}


}

