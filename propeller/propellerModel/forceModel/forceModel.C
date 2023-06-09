#include "forceModel.H"
#include "addToRunTimeSelectionTable.H"
#include "bladeGrid.H"
#include "linearInterpolation.H"


namespace Foam
{
    //set and define the type name "typeName"
    defineTypeNameAndDebug(forceModel,0);

    //Add to run time table to dynamically select the class based on
    //the dictionary propellerModel atribute
    addToRunTimeSelectionTable(propellerModel,forceModel,dictionary);


forceModel::forceModel
(
    const dictionary& dict
) : propellerModel(dict,typeName),
    gridDictionary_(dict.subDict("rotorGrid"))
{

    Info<<"Creating force Model"<<endl;
    rhoRef_ = dict.get<scalar>("rhoRef");
    control_ = fmControl::New(dict.subDict("control"),*this);

    List<scalar> J({0,1});
    List<scalar> CT({0,0.1});
    List<scalar> CQ({0,0.01});

    thrustCoeff_ = autoPtr<regularInterpolation<scalar,scalar,1>>
    ::NewFrom<linearInterpolation<scalar,scalar,1>>(FixedList<scalarList,1>({J}),CT);
    torqueCoeff_ = autoPtr<regularInterpolation<scalar,scalar,1>>
    ::NewFrom<linearInterpolation<scalar,scalar,1>>(FixedList<scalarList,1>({J}),CQ);


}

scalar forceModel::AxCoefficient(scalar thrust, scalar minRadius, scalar maxRadius)
{
    return 105.0/8.0*thrust/(constant::mathematical::pi*(3*minRadius+4*maxRadius)*(maxRadius-minRadius));
}

scalar forceModel::AthetaCoefficient(scalar torque, scalar minRadius, scalar maxRadius)
{
    return 105.0/8.0*torque/(constant::mathematical::pi*maxRadius*(3*minRadius+4*maxRadius)*(maxRadius-minRadius));
}

void forceModel::updateTensors()
{
    //Resize if num change
    if(gridTensor_.size()!=rotorGrid_->nCells())
    {
        gridTensor_.resize(rotorGrid_->nCells());
    }

    //Calculate all tensor
    forAll(gridTensor_,i)
    {
        gridTensor_[i]=cellBladeTensor(rotorGrid_->cells()[i]);
    }
}

tensor forceModel::cellBladeTensor(const gridCell &cell) const
{
    vector center = cell.center();
    return propellerModel::bladeTensor(rotorGrid_->geometry().cylindricalCS(),center,0,0);
}

scalar forceModel::getReferenceSpeed(const vectorField & U, const vector & normal)
{
    vector velAvg = average(U);
    return (-normal.inner(velAvg));
}

void forceModel::build(const rotorGeometry &rotorGeometry)
{
    rotorGrid_ = rotorGrid::New(gridDictionary_,rotorGeometry,*rotorFvMeshSel_,nullptr,1);
    this->updateTensors();
}

void forceModel::nextTimeStep(scalar dt)
{
    const bladeGrid* bg = dynamic_cast<const bladeGrid*>(rotorGrid_.get());
    if(bg)
    {
        //Update angle if its a bladeGrid
        scalar dpsi = dt*control_->getAngularVelocity();
        if(dpsi>constant::mathematical::twoPi)
        {
            Warning 
                << "Rotor angular step is bigger than 1 revolution: "
                << dpsi/constant::mathematical::twoPi
                << endl;
        }
        psi0_+= dpsi;
        if(psi0_>constant::mathematical::pi)
        {
            psi0_ -=constant::mathematical::twoPi;
        }

        List<scalar> initialPos = bg->getInitialAzimuth();
        forAll(initialPos,i)
        {
            initialPos[i] = initialPos[i]+psi0_;
        }
        rotorGrid_->setRotation(initialPos);
        this->updateTensors();
    }    
}

propellerResult forceModel::calculate(const vectorField &U, const scalarField *rhoField, volVectorField &force)
{
    propellerResult result;

    volVectorField fmForce(
        IOobject
        (
            "fmForce",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector(dimForce, Zero)
    );

    vector normal = rotorGrid_->geometry().direction();
    scalar speedRef = getReferenceSpeed(U,normal);

    scalar rhoRef = rhoField == nullptr ? rhoRef_ : average(*rhoField);
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
        // r-theta
        vector cellForce = cell.scaleForce(forceOverLen);

        // global
        vector globalForce = transform(localBlade,cellForce);

        //x-y-z local
        vector localForce = rotorGrid_->geometry().cartesianCS().localVector(globalForce);
        
        result.force += localForce;

        vector localPos = rotorGrid_->geometry().cartesianToCylindrical().globalPosition(cell.center());
        result.torque += localForce.cross(localPos);
        cell.applySource(force,cellVol,globalForce);
        cell.applyField<vector>(fmForce,cellForce);

    }

    scalar omega = control_->getAngularVelocity();

    result.update(omega, speedRef,rhoRef, maxR);

    if(rotorFvMeshSel_->mesh().time().writeTime())
    {
        fmForce.write();
    }

    return result;
}

propellerResult forceModel::calculate(const vectorField &U, const scalarField *rhoField) const
{
    propellerResult result;

    vector normal = rotorGrid_->geometry().direction();
    scalar speedRef = getReferenceSpeed(U,normal);

    scalar rhoRef = rhoField == nullptr ? rhoRef_ : average(*rhoField);
    scalar maxR = rotorGrid_->geometry().radius();
    scalar minR = rotorGrid_->geometry().innerRadius();
    scalar D = 2 * maxR;

    scalar J = control_->getJ();

    scalar Kt = thrustCoeff_->interpolate({J}).value();
    scalar Kq = torqueCoeff_->interpolate({J}).value();

    scalar thrust = Kt*rhoRef*pow(speedRef,2)*pow(D,2)/pow(J,2);
    scalar torque = Kq*rhoRef*pow(speedRef,2)*pow(D,3)/pow(J,2);

    scalar Ax = AxCoefficient(thrust,minR,maxR);
    scalar Atheta = AthetaCoefficient(torque,minR,maxR);

    const PtrList<gridCell>& cells = rotorGrid_->cells();

    forAll(cells,i)
    {
        const gridCell& cell = cells[i]; 
        const tensor& localBlade = gridTensor_[i];
        scalar cellR = cell.center().x();

        vector forceOverLen = ForceDistribution(Ax,Atheta,cellR,minR,maxR);
        // r-theta
        vector cellForce = cell.scaleForce(forceOverLen);

        // global
        vector globalForce = transform(localBlade,cellForce);

        //x-y-z local
        vector localForce = rotorGrid_->geometry().cartesianCS().localVector(globalForce);
        
        result.force += localForce;

        vector localPos = rotorGrid_->geometry().cartesianToCylindrical().globalPosition(cell.center());
        result.torque += localForce.cross(localPos);

    }

    scalar omega = control_->getAngularVelocity();

    result.update(omega, speedRef,rhoRef, maxR);

    return result;
}

vector forceModel::ForceDistribution(scalar Ax, scalar Atheta, scalar radius, scalar minRadius, scalar maxRadius)
{
    scalar r1h = minRadius/maxRadius;
    scalar r1=radius/maxRadius;
    scalar rStar = (r1-r1h)/(1-r1h);
    scalar sqrtRstar = sqrt(1-rStar);

    scalar fx = Ax*rStar*sqrtRstar;
    scalar ftheta = Atheta*rStar*sqrtRstar/(rStar*(1-r1h)+r1h);

    return vector(ftheta,0,fx);
}

}
 