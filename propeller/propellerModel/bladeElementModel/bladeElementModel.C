#include "bladeElementModel.H"
#include "addToRunTimeSelectionTable.H"
#include "rotorGrid.H"
#include "cubicSplineInterpolation.H"
#include "bladeGrid.H"
#include "bemControl.H"

namespace Foam
{
    //set and define the type name "typeName"
    defineTypeNameAndDebug(bladeElementModel,0);

    //Add to run time table to dynamically select the class based on
    //the dictionary propellerModel atribute
    addToRunTimeSelectionTable(propellerModel,bladeElementModel,dictionary);

const Enum<bladeElementModel::controlVar> 
bladeElementModel::controlVarNames_
({
    {controlVar::omega, "angularVelocity"},
    {controlVar::collectivePitch, "collectivePitch"},   
    {controlVar::ciclicPitchCos, "ciclicPitchCos"},   
    {controlVar::ciclicPitchSin, "ciclicPitchSin"}   
});

const Enum<bladeElementModel::outputVar>
bladeElementModel::outputVarNames_
({
    {outputVar::forceX, "forceX"},
    {outputVar::forceY, "forceY"},
    {outputVar::forceZ, "forceZ"},
    {outputVar::torqueX, "torqueX"},
    {outputVar::torqueY, "torqueY"},
    {outputVar::torqueZ, "torqueZ"},
    {outputVar::power, "power"},

});

bladeElementModel::bladeElementModel
(
    const dictionary& dict
) : 
    propellerModel(dict,typeName),
    airfoils_(dict.subDict("airfoils")),
    bladeModel_(airfoils_,dict.subDict("bladeModel")),
    gridDictionary(dict.subDict("rotorGrid")),

    polarCorrection_(dict)
{    
    dict.readEntry("nBlades",nBlades_);
    tipFactor_ = dict.getOrDefault<scalar>("tipFactor",1);

    rhoRef_ = dict.get<scalar>("rhoRef");
    speedRef_ = dict.get<scalar>("speedRef");
    nuRef_ = dict.get<scalar>("viscosity");
    soundSpeedRef_ = dict.get<scalar>("soundSpeed");

    control_ = bemControl::New(dict.subDict("control"),*this);
}


void bladeElementModel::build(const rotorGeometry& rotorGeometry)
{
    rotorGrid_ = rotorGrid::New(gridDictionary,rotorGeometry,*rotorFvMeshSel_,&bladeModel_,nBlades_);
    this->updateTensors();

    volScalarField selected
    (
        IOobject
        (
            "selectedCells",
            rotorFvMeshSel_->mesh().time().timeName(),
            rotorFvMeshSel_->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rotorFvMeshSel_->mesh(),
        dimensionedScalar(dimless, Zero)
    );
    volVectorField gridpos
    (
        IOobject
        (
            "gridPosition",
            rotorFvMeshSel_->mesh().time().timeName(),
            rotorFvMeshSel_->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rotorFvMeshSel_->mesh(),
        dimensionedVector(dimless, Zero)
    );

    forAll(rotorGrid_->cells(),i)
    {
        rotorGrid_->cells()[i].applyField<scalar>(selected.ref(false),1);
        rotorGrid_->cells()[i].applyField<vector>(gridpos.ref(false),rotorGrid_->cells()[i].center());
    }

    selected.write();
    gridpos.write();
    bladeModel_.writeBlade(300,"blade.csv");
}

propellerResult bladeElementModel::calculate(const vectorField& U, const scalarField* rhoField, volVectorField& force)
{
    propellerResult result;
    //Puntos de la discretizacion
    PtrList<gridCell>& cells = rotorGrid_->cells();
    
    scalar aoaMax = -VGREAT;
    scalar aoaMin = VGREAT;

    const scalarField& cellVol = rotorFvMeshSel_->mesh().V();

    volScalarField weights(
        IOobject
        (
            "weights",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar(dimless, Zero)
    );
    volVectorField bemForce(
        IOobject
        (
            "bemForce",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector(dimless, Zero)
    );
    volVectorField bladeLocalVel(
        IOobject
        (
            "bladeLocalVel",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector(dimVelocity, Zero)
    );
    volScalarField aoa(
        IOobject
        (
            "aoa",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar(dimless, Zero)
    );
    
    control_->correctControl(U,rhoField);

    scalar angularVelocity = control_->getAngularVelocity();
    //---CALCULATE VALUE ON INTEGRATION POINTS---//
    forAll(cells, i)
    {
        gridCell& cell = cells[i]; 
        const tensor& localBlade = gridTensor_[i];
        bemDebugData data;
        scalar rho = rhoField?(*rhoField)[i]:rhoRef_;
        vector forceOverLen = this->calculatePoint(U[i],rho,angularVelocity,cell,localBlade,data);
        vector cellforce = cell.scaleForce(forceOverLen);
        vector localForce = rotorGrid_->geometry().cartesianCS().localVector(cellforce);
        result.force += localForce;
        vector localPos = rotorGrid_->geometry().cartesianToCylindrical().globalPosition(cell.center());
        result.torque += localForce.cross(localPos);
        cell.applySource(force,cellVol,cellforce);
        cell.applyField<scalar>(weights,cell.weights());
        cell.applyField<vector>(bemForce,cellforce);
        cell.applyField<vector>(bladeLocalVel,transform(localBlade,control_->getBladeLocalVel(cell.azimuth(),cell.radius())));
        cell.applyField<scalar>(aoa,data.aoa);
        if(data.aoa>aoaMax) aoaMax=data.aoa;
        if(data.aoa<aoaMin) aoaMin=data.aoa;

    }

    result.update(angularVelocity,speedRef_,rhoRef_,rotorGrid_->geometry().radius());
    
    Info<< "- Max AoA: "<<aoaMax * 180/constant::mathematical::pi <<"º"<<endl;
    Info<< "- Min AoA: "<<aoaMin * 180/constant::mathematical::pi <<"º"<<endl;

    if(mesh().time().writeTime())
    {
        weights.write();
        bemForce.write();
        bladeLocalVel.write();
        aoa.write();
    }

    return result;
}

propellerResult bladeElementModel::calculate(const vectorField &U, const scalarField *rhoField) const
{
    propellerResult result;
    //Puntos de la discretizacion
    const PtrList<gridCell>& cells = rotorGrid_->cells();
    
    scalar angularVelocity = control_->getAngularVelocity();
    //---CALCULATE VALUE ON INTEGRATION POINTS---//
    forAll(cells, i)
    {
        const gridCell& cell = cells[i]; 
        const tensor& localBlade = gridTensor_[i];
        bemDebugData data;
        scalar rho = rhoField?(*rhoField)[i]:rhoRef_;
        vector forceOverLen = this->calculatePoint(U[i],rho,angularVelocity,cell,localBlade,data);
        vector cellforce = cell.scaleForce(forceOverLen);

        vector localForce = rotorGrid_->geometry().cartesianCS().localVector(cellforce);
        result.force += localForce;
        vector localPos = rotorGrid_->geometry().cartesianToCylindrical().globalPosition(cell.center());
        result.torque += localForce.cross(localPos);

    }

    result.update(angularVelocity,speedRef_,rhoRef_,rotorGrid_->geometry().radius());


    return result;
}

vector bladeElementModel::calculatePoint(const vector &U,scalar rho, scalar angularVelocity, const gridCell &cell, const tensor& bladeTensor, bemDebugData& data) const
{
    //---GET INTERPOLATED SECTION---//
    scalar radius = cell.radius();
    scalar azimuth = cell.azimuth();
    scalar maxRadius = rotorGrid_->geometry().radius();
    scalar chord,twist,sweep;
    interpolatedAirfoil airfoil;
    if(!bladeModel_.sectionAtRadius(radius/maxRadius,chord,twist,sweep,airfoil))
    {
        return Zero;
    }       

    //---GET RELATIVE VELOCITY---//
    //Local velocity vector
    vector localAirVel = invTransform(bladeTensor,U);
    
    //Get blade velocity
    vector relativeBladeVel = - control_->getBladeLocalVel(azimuth,radius);

    //Get relative air velocity
    vector relativeVel = localAirVel + relativeBladeVel;

    //y component is radial, thus "doesn't contribute to aerodinamic forces"
    relativeVel.y()=0;
    scalar relativeSpeed = mag(relativeVel);

    //Airspeed angle (positive when speed is from "below" airfoil)
    scalar phi = bladeElementModel::AngleOfIncidence(relativeVel);

    //Get collective and ciclic pitch
    scalar pitch = control_->getPitch(azimuth);

    //Angle of atack
    scalar AoA = twist + pitch - phi;

    //---COMPUTE AERODYNAMIC FORCE ON BLADE---//
    scalar re = rho*relativeSpeed*chord/nuRef_;
    scalar mach = relativeSpeed/soundSpeedRef_;

    scalar cl = airfoil.cl(AoA,re,mach);
    scalar cd = airfoil.cd(AoA,re,mach);
    polarCorrection_.correct(cl,cl,AoA,chord,radius,maxRadius,angularVelocity,mag(localAirVel));
    //Add tip factor effect:
    if(radius/maxRadius>=tipFactor_)
    {
        cl=0.0;
    }
    
    //Calculate aerodinamic forces
    scalar dynamicPressure = 0.5 * rho * chord * relativeSpeed * relativeSpeed;
    scalar lift = cl * dynamicPressure;
    scalar drag = cd * dynamicPressure;

    //Open foam code is not decomposing Lift and drag
    //from wind axis to blade axis
    //Project over normal components
    vector localForce(lift*sin(phi) + drag * cos(phi),0,lift * cos(phi) - drag * sin(phi));
    
    
    //Back to global ref frame
    vector globalForceOverLenght = transform(bladeTensor,localForce);

    data.cl = cl;
    data.cd = cd;
    data.phi = phi;
    data.aoa = AoA;
    data.radius = radius;
    data.chord = chord;
    data.localVel = U;
    return globalForceOverLenght;
}
bool bladeElementModel::nextTimeStep(scalar dt)
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
            initialPos[i] = (control_->getAzimuth(initialPos[i]+psi0_));
        }
        bg->setRotation(initialPos);
        this->updateTensors();

        return true;
    }
    return false;
}
inline scalar bladeElementModel::AngleOfIncidenceSTAR(const vector &relativeLocalVel)
{
    return atan2(-relativeLocalVel.z(),sign(relativeLocalVel.x())*sqrt(pow(relativeLocalVel.x(),2)+pow(relativeLocalVel.y(),2))); //as starccm+
}
inline scalar bladeElementModel::AngleOfIncidence(const vector &relativeLocalVel)
{
    return atan2(-relativeLocalVel.z(),relativeLocalVel.x());
}

bool bladeElementModel::isTimeAcuratte() const
{
    return dynamic_cast<const bladeGrid*>(rotorGrid_.get()) != nullptr;
}

void bladeElementModel::updateTensors()
{
    volVectorField x
    (
        IOobject
        (
            "xAxis",
            rotorFvMeshSel_->mesh().time().timeName(),
            rotorFvMeshSel_->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rotorFvMeshSel_->mesh(),
        dimensionedVector(dimless, Zero)
    );
    volVectorField y
    (
        IOobject
        (
            "yAxis",
            rotorFvMeshSel_->mesh().time().timeName(),
            rotorFvMeshSel_->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rotorFvMeshSel_->mesh(),
        dimensionedVector(dimless, Zero)
    );
    volVectorField z
    (
        IOobject
        (
            "zAxis",
            rotorFvMeshSel_->mesh().time().timeName(),
            rotorFvMeshSel_->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rotorFvMeshSel_->mesh(),
        dimensionedVector(dimless, Zero)
    );

    //Resize if num change
    if(gridTensor_.size()!=rotorGrid_->nCells())
    {
        gridTensor_.resize(rotorGrid_->nCells());
    }

    forAll(gridTensor_,i)
    {
        gridTensor_[i]=cellBladeTensor(rotorGrid_->cells()[i]);
        rotorGrid_->cells()[i].applyField<vector>(x,gridTensor_[i].col<0>());
        rotorGrid_->cells()[i].applyField<vector>(y,gridTensor_[i].col<1>());
        rotorGrid_->cells()[i].applyField<vector>(z,gridTensor_[i].col<2>());
    }
    x.write();
    y.write();
    z.write();

}

tensor bladeElementModel::cellBladeTensor(const gridCell &cell) const
{
    vector center = cell.center();
    scalar radius = center.x();
    scalar azimuth = center.y();
    scalar sweep = bladeModel_.getSweep(radius);
    scalar flapping = control_->getFlapping(azimuth);

    return propellerModel::bladeTensor(rotorGrid_->geometry().cylindricalCS(),center,flapping,sweep);
}

}
