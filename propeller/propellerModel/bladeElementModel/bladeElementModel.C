#include "bladeElementModel.H"
#include "bladeSection.H"
#include "addToRunTimeSelectionTable.H"


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
    airfoils_(dict.subDict("airfoils")),
    bladeModel_(airfoils_,dict.subDict("bladeModel")),
    propellerModel(dict,typeName)
{
    Info<<"Creating blade Element Model"<<endl;
    
}

Foam::scalar Foam::bladeElementModel::radius() const
{
    return bladeModel_.maxRadius();
}

void Foam::bladeElementModel::build(const rotorGeometry& rotorGeometry)
{
    rotorDiscrete_.buildCoordinateSystem(rotorGeometry);
    rotorDiscrete_.fromRotorMesh(*rotorMesh_);

    bladeModel_.setMaxRadius(rotorGeometry.radius);
}


Foam::propellerResult Foam::bladeElementModel::calculate(const vectorField& U,volVectorField& force)
{
    propellerResult result;
    //Puntos de la discretizacion
    const List<vector> cylPoints = rotorDiscrete_.cylPoints();
    //Tensor de cada punto local to global
    const List<tensor> bladeCS = rotorDiscrete_.localBladeCS();

    //Velocidad angular
    double rpm = 1000;
    double pi = Foam::constant::mathematical::pi;
    double omega = rpm*pi/30;
    
    List<scalar> aoaList(cylPoints.size());

    volScalarField aoaField
    (
        IOobject
        (
            "propeller:AoA",
            rotorMesh_->mesh().time().timeName(),
            rotorMesh_->mesh()
        ),
        rotorMesh_->mesh(),
        dimensionedScalar(dimless, Zero)
    );

    forAll(cylPoints, i)
    {
        //Get local radius
        scalar radius = cylPoints[i].x();
        auto bladeSec = bladeModel_.sectionAtRadius(radius);
        scalar chord = bladeSec.chord();

        if(chord == 0)
        {
             continue;
        }
        scalar twist = bladeSec.twist();
        scalar n_blade = 2;
        scalar average_fact = n_blade / (2 * pi * radius);

        //Get cell area and volum
        scalar area = rotorMesh_->areas()[i];
        label celli = rotorMesh_->cells()[i];
        scalar volume = rotorMesh_->mesh().V()[celli];
        
        
        //Local rotation tensor
        const tensor& bladeTensor = bladeCS[i];

        //Global coordinate vector
        vector airVel = U[i];
        vector localAirVel = invTransform(bladeTensor,airVel);
        

        //Get blade velocity
        vector relativeBladeVel(omega*radius,0,0);

        //Get relative air velocity
        vector relativeVel = localAirVel + relativeBladeVel;

        //y component is radial, thus "doesn't contribute to aerodinamic forces"
        relativeVel.y()=0;
        scalar relativeSpeed = mag(relativeVel);

        //Info<<"Rel local vel: "<<relativeSpeed<<endl;
        //Airspeed angle (positive when speed is from "below" airfoil)
        scalar phi = atan2(-relativeVel.z(),relativeVel.x());
        //Info<<"Phi: "<<phi<<endl;

        

        //Angle of atack
        scalar AoA = twist - phi;

        scalar rho = this->refRho;
        scalar nu = 1e-5;
        scalar re = rho*relativeSpeed*chord/nu;
        scalar c = 345;
        scalar mach = relativeSpeed/c;

        scalar cl = bladeSec.cl(AoA,re,mach);
        scalar cd = bladeSec.cd(AoA,re,mach);

        aoaField[celli] = AoA;
        aoaList[i]=AoA;
       

        //Calculate aerodinamic forces
        scalar lift = average_fact * 0.5 * rho * cl * chord * relativeSpeed * relativeSpeed * area;
        scalar drag = average_fact * 0.5 * rho * cd * chord * relativeSpeed * relativeSpeed * area;
        
        
        //Project over normal components
        vector normalForce(0,0,lift * cos(phi) - drag * sin(phi));
        vector tangentialForce(lift*sin(phi) + drag * cos(phi),0,0);
        
        vector totalAerForce = normalForce + tangentialForce;

        result.force += totalAerForce;
        result.torque += (tangentialForce.x() * radius);

        //Back to global ref frame
        totalAerForce = transform(bladeTensor,totalAerForce);

        //Add source term
        force[celli] = -totalAerForce / volume;
        
    }
    if(rotorMesh_->mesh().time().writeTime())
    {
        aoaField.write();
    }

    //TODO: obtain real time step
    //rotorDynamics_.integrate(-totalMoment,0.005);
    result.power = result.torque * omega;
    result.eta=result.force.z() * this->refV / result.power;
    result.J = this->refV/(omega*rotorDiscrete_.geometry().radius);

    /*Info<< "Total Lift: "<<totalLift<<endl;
    Info<< "Total Drag: "<<totalDrag<<endl;
    Info<< "Total thrust: "<<totalThrust<<endl;
    Info<< "Total power:" <<power <<endl;
    Info<< "Total Moment z: " <<totalMoment<<endl;
    Info<< "RPM: "<< omega * 30/pi<<endl;
    Info<< "Rad/s: " <<omega<<endl;
    Info<< "Max AoA: "<< max(aoaList) * 180/pi<<endl;
    Info<< "Min AoA: "<< min(aoaList) * 180/pi<<endl;
    Info<< "Max Sampled Vel: "<<max(U)<<endl;*/
    return result;

}
