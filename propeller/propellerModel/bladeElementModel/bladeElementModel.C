#include "bladeElementModel.H"
#include "addToRunTimeSelectionTable.H"
#include "rotorGrid.H"
#include "cubicSplineInterpolation.H"
#include "bladeGrid.H"
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
    control_(dict.subDict("control"))
{    
    dict.readEntry("nBlades",nBlades_);
    tipFactor_ = dict.getOrDefault<scalar>("tipFactor",1);
}


void Foam::bladeElementModel::build(const rotorGeometry& rotorGeometry)
{
    rotorGrid_ = rotorGrid::New(gridDictionary,rotorGeometry,*rotorFvMeshSel_,bladeModel_,nBlades_);

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

Foam::propellerResult Foam::bladeElementModel::calculate(const vectorField& U,scalar angularVelocity, volVectorField& force, scalar theta)
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
    Info<<"theta: "<<theta<<endl;
    const bladeGrid* bg = dynamic_cast<const bladeGrid*>(rotorGrid_.get());
    if(bg)
    {
        List<scalar> initialPos = bg->getInitialAzimuth();
        forAll(initialPos,i)
        {
            initialPos[i] = (control_.getAzimuth(initialPos[i]+theta));
        }
        rotorGrid_->setRotation(initialPos);
    }


    //---CALCULATE VALUE ON INTEGRATION POINTS---//
    forAll(cells, i)
    {
        gridCell& cell = cells[i]; 
        bemDebugData data;
  
        vector forceOverLen = this->calculatePoint(U[i],angularVelocity,cell,data);


        vector cellforce = cell.scaleForce(forceOverLen);

        vector localForce = rotorGrid_->geometry().cartesianCS().localVector(cellforce);
        result.force += localForce;
        vector localPos = rotorGrid_->geometry().cartesianToCylindrical().globalPosition(cell.center());
        result.torque += localForce.cross(localPos);
        cell.applySource(force,cellVol,cellforce);
        cell.applyField<scalar>(weights,cell.weights());
        cell.applyField<vector>(bemForce,cellforce);
        cell.applyField<vector>(bladeLocalVel,transform(cell.localBlade(),control_.getBladeLocalVel(cell.azimuth(),angularVelocity,cell.radius())));
        cell.applyField<scalar>(aoa,data.aoa);
        if(data.aoa>aoaMax) aoaMax=data.aoa;
        if(data.aoa<aoaMin) aoaMin=data.aoa;

    }

    reduce(aoaMax,maxOp<scalar>());
    reduce(aoaMin,minOp<scalar>());

    reduce(result.force,sumOp<vector>());
    reduce(result.torque,sumOp<vector>());

    result.power = result.torque.z() * angularVelocity;

    result.updateEta(this->refV);
    //Use J definition Vref/(rps * D)
    result.updateJ(this->refV,angularVelocity,rotorGrid_->geometry().radius().get());
    result.updateCT(this->refRho,angularVelocity,rotorGrid_->geometry().radius().get());
    result.updateCP(this->refRho,angularVelocity,rotorGrid_->geometry().radius().get());
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

Foam::vector Foam::bladeElementModel::calculatePoint(const vector &U, scalar angularVelocity, const gridCell &cell, bemDebugData& data)
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
    //Local rotation tensor
    tensor bladeTensor = cell.localBlade();

    //---GET RELATIVE VELOCITY---//
    //Local velocity vector
    vector localAirVel = invTransform(bladeTensor,U);
    
    //Get blade velocity
    vector relativeBladeVel = - control_.getBladeLocalVel(azimuth,angularVelocity,radius);

    //Get relative air velocity
    vector relativeVel = localAirVel + relativeBladeVel;

    //y component is radial, thus "doesn't contribute to aerodinamic forces"
    //relativeVel.y()=0;
    scalar relativeSpeed = mag(relativeVel);

    //Airspeed angle (positive when speed is from "below" airfoil)
    scalar phi = bladeElementModel::AngleOfIncidenceSTAR(relativeVel);

    //Get collective and ciclic pitch
    scalar pitch = control_.getPitch(azimuth);

    //Angle of atack
    scalar AoA = twist + pitch - phi;

    //---COMPUTE AERODYNAMIC FORCE ON BLADE---//
    scalar rho = this->refRho;
    scalar nu = 1e-5;
    scalar re = rho*relativeSpeed*chord/nu;
    scalar c = 345;
    scalar mach = relativeSpeed/c;

    scalar cl = airfoil.cl(AoA,re,mach);
    scalar cd = airfoil.cd(AoA,re,mach);

    //Add tip factor effect:
    if(radius/rotorGrid_->geometry().radius().get()>=tipFactor_)
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
inline Foam::scalar Foam::bladeElementModel::AngleOfIncidenceSTAR(const vector &relativeLocalVel)
{
    return atan2(-relativeLocalVel.z(),sign(relativeLocalVel.x())*sqrt(pow(relativeLocalVel.x(),2)+pow(relativeLocalVel.y(),2))); //as starccm+
    //return atan2(-relativeVel.z(),relativeVel.x());
}
inline Foam::scalar Foam::bladeElementModel::AngleOfIncidence(const vector &relativeLocalVel)
{
    //return atan2(-relativeLocalVel.z(),sign(relativeLocalVel.x())*sqrt(pow(relativeLocalVel.x(),2)+pow(relativeLocalVel.y(),2))); //as starccm+
    return atan2(-relativeLocalVel.z(),relativeLocalVel.x());
}
