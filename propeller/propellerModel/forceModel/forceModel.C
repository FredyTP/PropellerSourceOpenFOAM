#include "forceModel.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{
    //set and define the type name "typeName"
    defineTypeNameAndDebug(forceModel,0);

    //Add to run time table to dynamically select the class based on
    //the dictionary propellerModel atribute
    //addToRunTimeSelectionTable(propellerModel,forceModel,dictionary);
}

Foam::forceModel::forceModel
(
    const dictionary& dict
) : propellerModel(dict,typeName),
    gridDictionary_(dict.subDict("rotorGrid")),
{

    Info<<"Creating force Model"<<endl;
    rhoRef_ = dict.get<scalar>("rhoRef");

    control_ = fmControl::New(dict.subDict("control"),*this);
}

scalar Foam::forceModel::AxCoefficient(scalar thrust, scalar minRadius, scalar maxRadius)
{
    return 105.0/8.0*thrust/(constant::mathematical::pi*(3*minRadius+4*maxRadius)*(maxRadius-minRadius));
}

scalar Foam::forceModel::AthetaCoefficient(scalar torque, scalar minRadius, scalar maxRadius)
{
    return 105.0/8.0*torque/(constant::mathematical::pi*maxRadius*(3*minRadius+4*maxRadius)*(maxRadius-minRadius));
}

propellerResult Foam::forceModel::calculate(const vectorField &U, const scalarField *rhoField, volVectorField &force)
{
    propellerResult result;

    vector normal = rotorGrid_->geometry().direction();
    scalar speedRef = normal.inner(cmptAv(U));
    scalar rhoRef = rhoField == nullptr ? rhoRef_ : cmptAv(*rhoField);
    scalar maxR = rotorGrid_->geometry().radius();
    scalar minR = rotorGrid_->geometry().innerRadius();
    scalar D = 2 * maxR;

    control_->correctControl(U,rhoField);
    scalar J = control_->getJ();

    scalar Kt = thrustCoeff_->interpolate({J}).value();
    scalar Kq = torqueCoeff_->interpolate({J}).value();

    scalar thrust = Kt*rhoRef*pow(speedRef,2)*pow(D,2)/pow(J,2);
    scalar torque = Kq*rhoRef*pow(speedRef,2)*pow(D,3)/pow(J,2);

    Info<<"Table thrust: "<<thrust<<endl;
    Info<<"Table torque: "<<torque<<endl;

    scalar Ax = AxCoefficient(thrust,minR,maxR);
    scalar Atheta = AthetaCoefficient(torque,minR,maxR);

    PtrList<gridCell>& cells = rotorGrid_->cells();
    const scalarField& cellVol = rotorFvMeshSel_->mesh().V();

    forAll(cells,i)
    {
        gridCell& cell = cells[i]; 
        const tensor& localBlade = gridTensor_[i];
        scalar cellR = cell.center().x();

        vector forceOverLen = ForceDistribution(Ax,Atheta,cellR,minR,maxR);
        vector cellForce = cell.scaleForce(forceOverLen);
        vector globalForce = transform(localBlade,localForce);
        vector localForce = rotorGrid_->geometry().cartesianCS().localVector(globalForce);
        result.force += localForce;
        vector localPos = rotorGrid_->geometry().cartesianToCylindrical().globalPosition(cell.center());
        result.torque += localForce.cross(localPos);
        cell.applySource(force,cellVol,globalForce);
        cell.applyField<vector>(fmForce,cellforce);

    }

    scalar omega = control_->getAngularVelocity();
    result.update(omega, speedRef,rhoRef, maxR);

    return result;
}

vector Foam::forceModel::ForceDistribution(scalar Ax, scalar Atheta, scalar radius, scalar minRadius, scalar maxRadius)
{
    scalar r1h = minRadius/maxRadius;
    scalar r1=radius/maxRadius;
    scalar rStar = (r1-r1h)/(1-r1h);
    scalar sqrtRstar = sqrt(1-rStar);

    scalar fx = Ax*rStar*sqrtRstar
    scalar ftheta = Atheta*rStar*sqrtRstar/(rStar*(1-r1h)+r1h);

    return vector(fx,ftheta,0);

}
 