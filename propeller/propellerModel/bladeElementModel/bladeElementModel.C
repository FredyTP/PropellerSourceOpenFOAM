#include "bladeElementModel.H"
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
}

void Foam::bladeElementModel::calculate(volVectorField& force)
{
    scalar totalLift = 0;
    scalar totalDrag = 0;
    scalar totalThrust = 0;

    //Puntos de la discretizacion
    const List<vector> cylPoints = rotorDiscrete_.cylPoints();
    //Tensor de cada punto local to global
    const List<tensor> bladeCS = rotorDiscrete_.localBladeCS();

    //Velocidad angular
    double rpm = 1000;
    double omega = rpm*Foam::constant::mathematical::twoPi/60;
    double pi = Foam::constant::mathematical::pi;

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
        scalar chord = bladeModel_.chordAtRadius(radius);
        scalar twist = bladeModel_.twistAtRadius(radius);
        scalar n_blade = 2;
        scalar average_fact = n_blade / (2 * pi * radius);

        const tensor& bladeTensor = bladeCS[i];

        //Global coordinate vector
        vector airVel(0,0,0);
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

        //Get cell area and volum
        scalar area = rotorMesh_->areas()[i];
        label celli = rotorMesh_->cells()[i];
        scalar volume = rotorMesh_->mesh().V()[celli];

        //Angle of atack
        scalar AoA = twist - phi;
       
        scalar rho = 1.225;
        scalar nu = 1e-5;
        scalar re = rho*relativeSpeed*chord/nu;
        scalar c = 345;
        scalar mach = relativeSpeed/c;

        scalar cl = airfoils_.getAirfoil(0)->cl(AoA,re,mach);
        scalar cd = airfoils_.getAirfoil(0)->cd(AoA,re,mach);

        aoaField[celli] = AoA;
       

        //Calculate aerodinamic forces
        scalar lift = average_fact * 0.5 * rho * cl * chord * relativeSpeed * relativeSpeed * area;
        scalar drag = average_fact * 0.5 * rho * cd * chord * relativeSpeed * relativeSpeed * area;
        
        totalLift += lift;
        totalDrag += drag;
    
        //Project over normal components
        vector normalForce(0,0,lift * cos(phi) - drag * sin(phi));
        vector tangentialForce(lift*sin(phi) + drag * cos(phi),0,0);
        totalThrust += normalForce.z();
        //Back to global ref frame
        vector totalAerForce = normalForce + tangentialForce;
        totalAerForce = transform(bladeTensor,totalAerForce);

        //Add source term
        force[celli] = -totalAerForce / volume;
        
    }
    if(rotorMesh_->mesh().time().writeTime())
    {
        aoaField.write();
    }

    Info<< "Total Lift: "<<totalLift<<endl;
    Info<< "Total Drag: "<<totalDrag<<endl;
    Info<< "Total thrust: "<<totalThrust<<endl;

}