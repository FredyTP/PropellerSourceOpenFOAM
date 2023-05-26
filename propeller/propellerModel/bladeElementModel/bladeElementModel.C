#include "bladeElementModel.H"
#include "addToRunTimeSelectionTable.H"
#include "rotorGrid.H"
#include "cubicSplineInterpolation.H"
#include "bladeGrid.H"
#include "axisAngleRotation.H"
namespace Foam
{
    //set and define the type name "typeName"
    defineTypeNameAndDebug(bladeElementModel,0);

    //Add to run time table to dynamically select the class based on
    //the dictionary propellerModel atribute
    addToRunTimeSelectionTable(propellerModel,bladeElementModel,dictionary);
}

Foam::bladeElementModel::bladeElementModel
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


void Foam::bladeElementModel::build(const rotorGeometry& rotorGeometry)
{
    rotorGrid_ = rotorGrid::New(gridDictionary,rotorGeometry,*rotorFvMeshSel_,bladeModel_,nBlades_);
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

Foam::propellerResult Foam::bladeElementModel::calculate(const vectorField& U, const scalarField* rhoField, volVectorField& force)
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

    reduce(aoaMax,maxOp<scalar>());
    reduce(aoaMin,minOp<scalar>());

    reduce(result.force,sumOp<vector>());
    reduce(result.torque,sumOp<vector>());

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

Foam::propellerResult Foam::bladeElementModel::calculate(const vectorField &U, const scalarField *rhoField) const
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

    reduce(result.force,sumOp<vector>());
    reduce(result.torque,sumOp<vector>());

    result.update(angularVelocity,speedRef_,rhoRef_,rotorGrid_->geometry().radius());


    return result;
}

Foam::vector Foam::bladeElementModel::calculatePoint(const vector &U,scalar rho, scalar angularVelocity, const gridCell &cell, const tensor& bladeTensor, bemDebugData& data) const
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
void Foam::bladeElementModel::nextTimeStep(scalar dt)
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
            initialPos[i] = (control_->getAzimuth(initialPos[i]+psi0_));
        }
        rotorGrid_->setRotation(initialPos);
        this->updateTensors();
    }
}
inline Foam::scalar Foam::bladeElementModel::AngleOfIncidenceSTAR(const vector &relativeLocalVel)
{
    return atan2(-relativeLocalVel.z(),sign(relativeLocalVel.x())*sqrt(pow(relativeLocalVel.x(),2)+pow(relativeLocalVel.y(),2))); //as starccm+
}
inline Foam::scalar Foam::bladeElementModel::AngleOfIncidence(const vector &relativeLocalVel)
{
    return atan2(-relativeLocalVel.z(),relativeLocalVel.x());
}

bool Foam::bladeElementModel::isTimeAcuratte() const
{
    return dynamic_cast<const bladeGrid*>(rotorGrid_.get()) != nullptr;
}

void Foam::bladeElementModel::updateTensors()
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

    gridTensor_.resize(rotorGrid_->nCells());
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

Foam::tensor Foam::bladeElementModel::cellBladeTensor(const gridCell &cell) const
{
    vector center = cell.center();
    scalar radius = center.x();
    scalar azimuth = center.y();
    scalar sweep = bladeModel_.getSweep(radius);
    scalar flapping = control_->getFlapping(azimuth);

    return bladeTensor(rotorGrid_->geometry().cylindricalCS(),center,flapping,sweep);
}

Foam::tensor Foam::bladeElementModel::bladeTensor(const coordSystem::cylindrical &cylCS, const point &localPoint, scalar flapping, scalar sweep)
{
    // z- up, y -outwards from center, x perpendicular y,z (leading edge to trailing edge)
    point global, origin;

    global = cylCS.globalPosition(localPoint);
    origin = cylCS.origin();

    tensor rotTensor(cylCS.R());

    // z-axis
    const vector ax3 = rotTensor.col<2>(); // == e3 (already normalized)

    // y-axis (radial direction)
    vector ax2(global - origin);

    ax2.removeCollinear(ax3);

    const scalar magAxis2(mag(ax2));

    // Trap zero size and colinearity
    if (magAxis2 < SMALL)
    {
        return rotTensor;
    }

    ax2 /= magAxis2; // normalise

    // Replace with updated local axes

    rotTensor.col<0>(ax2 ^ ax3);
    rotTensor.col<1>(ax2);
    //Rotate around X axis to get the flapping coordsys
    tensor rotX = coordinateRotations::axisAngle::rotation(vector::components::X,flapping,false);
    //Rotate around new Z axis to get the sweeped coordSys
    tensor rotZ = coordinateRotations::axisAngle::rotation(vector::components::Z,sweep,false);
    Info<<"flapping"<<flapping<<endl;
    return rotTensor.inner(rotX);
}
