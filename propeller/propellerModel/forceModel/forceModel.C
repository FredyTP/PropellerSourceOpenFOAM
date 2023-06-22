#include "forceModel.H"
#include "addToRunTimeSelectionTable.H"
#include "bladeGrid.H"
#include "LinearInterpolation.H"
#include "fmControl.H"
#include "readHelper.H"

namespace Foam
{
    //set and define the type name "typeName"
    defineTypeNameAndDebug(forceModel,0);

    //Add to run time table to dynamically select the class based on
    //the dictionary propellerModel atribute
    addToRunTimeSelectionTable(propellerModel,forceModel,dictionary);

const Enum<forceModel::controlVar> 
forceModel::controlVarNames_
({
    {controlVar::omega, "angularVelocity"}
});

const Enum<forceModel::outputVar>
forceModel::outputVarNames_
({
    {outputVar::forceZ, "forceZ"},
    {outputVar::torqueZ, "torqueZ"},
    {outputVar::power, "power"}
});

forceModel::forceModel
(
    word sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)  : 
    propellerModel(sourceName,dict,mesh),
    gridDictionary_(dict.subDict("rotorGrid"))
{
    read(dict);
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
scalar forceModel::getReferenceSpeed(const vectorField & U) const
{
    vector normal = rotorGrid_->geometry().direction();
    return forceModel::getReferenceSpeed(U,normal);
}
scalar forceModel::getReferenceSpeed(const vectorField & U, const vector & normal)
{
    vector velAvg = average(U);
    return (-normal.inner(velAvg));
}

void forceModel::read(const dictionary &dict)
{
    Info<<"Creating force Model"<<endl;
    rhoRef_ = dict.get<scalar>("rhoRef");
    control_ = fmControl::New(dict.subDict("control"),*this);

    const dictionary& curvesDict = dict.subDict("curves");
    thrustCoeff_ = util::NewInterpolationFromDict(curvesDict,"J","CT");
    torqueCoeff_ = util::NewInterpolationFromDict(curvesDict,"J","CQ");
}

void forceModel::build(const rotorGeometry &rotorGeometry)
{
    rotorGrid_ = rotorGrid::New(gridDictionary_,rotorGeometry,*rotorFvMeshSel_,nullptr,1);
    this->updateTensors();
}

bool forceModel::nextTimeStep(scalar dt)
{
    bladeGrid* bg = dynamic_cast<bladeGrid*>(rotorGrid_.get());
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
        bg->setRotation(initialPos);
        this->updateTensors();
        return true;
    }    
    return false;
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


    scalar speedRef = getReferenceSpeed(U);

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
        scalar cellR = cell.getCellCenter().x();

        scalar dr2 = cell.dr()/2;
        vector forceOverLen = ForceIntergral(Ax,Atheta,cellR-dr2,cellR+dr2,minR,maxR)/(cell.dr());
        // r-theta
        vector cellForce = cell.scaleForce(forceOverLen);

        // global
        vector globalForce = transform(localBlade,cellForce);

        //x-y-z local
        vector localForce = rotorGrid_->geometry().cartesianCS().localVector(globalForce);
        
        result.force += localForce;
        vector adimCenter = cell.center();
        adimCenter.x()=(adimCenter.x()-minR)/(maxR-minR);
        vector localPos = rotorGrid_->geometry().cartesianToCylindrical().globalPosition(adimCenter);
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
        scalar cellR = cell.getCellCenter().x();
        scalar dr2 = cell.dr()/2;
        vector forceOverLen = ForceIntergral(Ax,Atheta,cellR-dr2,cellR+dr2,minR,maxR)/(cell.dr());
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

vector forceModel::ForceIntergralFunction(scalar Ax, scalar Atheta, scalar radius, scalar minRadius, scalar maxRadius)
{
    scalar r1h = minRadius/maxRadius;
    scalar r1=radius/maxRadius;
    scalar rStar = (r1-r1h)/(1-r1h);
    
    //Force X
    scalar fx = -2*Ax*(3*rStar+2)*pow(1-rStar,1.5)/15;
    
    //Force theta
    scalar sqrt1h=sqrt(1-r1h);
    scalar sqrtr = sqrt(1-rStar);

    scalar a = (2.0*sqrtr*(1-rStar+r1h*(2+rStar))/(3.0*pow(r1h-1,2)));
    scalar den = pow(sqrt1h,5);

    scalar b = r1h*log((1-sqrt1h*sqrtr)/(1+sqrt1h*sqrtr))/den;
    scalar ftheta = -Atheta*(a+b);

    return vector(ftheta,0,fx);
}
scalar forceModel::AxCoefficient(scalar thrust, scalar minRadius, scalar maxRadius)
{
    return thrust*15.0/4.0;
}

scalar forceModel::AthetaCoefficient(scalar torque, scalar minRadius, scalar maxRadius)
{
    scalar h = minRadius/maxRadius;
    if(h==0.0)
    {
        return torque*3/2;
    }
    scalar a = (2.0*(-2+9*h+8*pow(h,2))/(15.0*pow(h-1,3)));
    scalar sqrt1h=sqrt(1-h);
    scalar den = pow(sqrt1h,7);

    scalar b = pow(h,2)*log((1-sqrt1h)/(1+sqrt1h))/den;

    return torque/(a-b);
}
vector forceModel::ForceIntergral(scalar Ax, scalar Atheta, scalar radius0, scalar radius1, scalar minRadius, scalar maxRadius)
{
    //radius1 > radius0
    return 
        ForceIntergralFunction(Ax,Atheta,radius1,minRadius,maxRadius)-
        ForceIntergralFunction(Ax,Atheta,radius0,minRadius,maxRadius);
}

}
 