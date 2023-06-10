#include "fmTargetValue.H"
#include "addToRunTimeSelectionTable.H"
#include "forceModel.H"
#include "functionSolver.H"
namespace Foam
{
 
defineTypeNameAndDebug(fmTargetValue,0);

addToRunTimeSelectionTable(fmControl,fmTargetValue,dictionary);

fmTargetValue::fmTargetValue(const dictionary& dict,const forceModel& fmModel)
 : fmControl(dict,fmModel)
{
    targetThrust_ = dict.readIfPresent("thrust",thrust_);
    if(!targetThrust_)
    {
        torque_ = dict.readEntry("torque",torque_);
    }
}

void fmTargetValue::correctControl(const vectorField & U, const scalarField * rhoField)
{
    vector normal = forceModel_.grid()->geometry().direction();
    scalar diameter = 2 * forceModel_.grid()->geometry().radius();
    scalar speedRef = forceModel::getReferenceSpeed(U,normal);
    scalarField x0(1,0.5);

    if(targetThrust_)
    {
        auto trim = util::functionSolver::NewtonRapson(1,forceFunction(U,rhoField),x0,x0/100,0.5,100,1e-12,true);

        Info<< "Final force: "<<forceFunction(U,rhoField)(scalarList(1,J_))+thrust_<<endl;
    }
    else
    {

    }
}
scalar fmTargetValue::getAngularVelocity() const
{
    return omega_;
}
scalar fmTargetValue::getJ() const
{
    return J_;
}

std::function<scalarField(scalarField)> fmTargetValue::forceFunction(const vectorField &U, const scalarField *rhoField)
{
    return [this,U,rhoField](const scalarField& J)
    {
        this->J_ = J[0];
        vector normal = forceModel_.grid()->geometry().direction();
        scalar speedRef = forceModel::getReferenceSpeed(U,normal);
        this->omega_ = this->getOmegaFromJ(this->J_,speedRef);

        auto result = this->forceModel_.calculate(U,rhoField);

        return scalarField(1,result.thrust()-this->thrust_); 
    };
}
}

