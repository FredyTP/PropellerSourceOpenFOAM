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

    bladeModel_.setMaxRadius(rotorGeometry.radius());
}


Foam::propellerResult Foam::bladeElementModel::calculate(const vectorField& U,scalar angularVelocity, volVectorField& force)
{
    propellerResult result;
    //Puntos de la discretizacion
    const List<vector> cylPoints = rotorDiscrete_.cylPoints();
    //Tensor de cada punto local to global
    const List<tensor> bladeCS = rotorDiscrete_.localBladeCS();

    //Velocidad angular
    scalar pi = Foam::constant::mathematical::pi;
    const scalar omega = angularVelocity;
    
    List<scalar> aoaList(cylPoints.size());
    List<vector> pressOnPoints(cylPoints.size(),vector(0,0,0));

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

    volScalarField radiusField
    (
        IOobject
        (
            "propeller:radius",
            rotorMesh_->mesh().time().timeName(),
            rotorMesh_->mesh()
        ),
        rotorMesh_->mesh(),
        dimensionedScalar(dimless, Zero)
    );

    volScalarField clField
    (
        IOobject
        (
            "propeller:cl",
            rotorMesh_->mesh().time().timeName(),
            rotorMesh_->mesh()
        ),
        rotorMesh_->mesh(),
        dimensionedScalar(dimless, Zero)
    );


    volScalarField cdField
    (
        IOobject
        (
            "propeller:cd",
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
        //scalar area = rotorMesh_->areas()[i];
        //label celli = rotorMesh_->cells()[i];
        //scalar volume = rotorMesh_->mesh().V()[celli];
        
        
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

        //clField[celli]=cl;
        //cdField[celli]=cd;
        //aoaField[celli] = AoA;
        //aoaList[i]=AoA;
       

        //Calculate aerodinamic forces
        scalar lift = average_fact * 0.5 * rho * cl * chord * relativeSpeed * relativeSpeed;
        scalar drag = average_fact * 0.5 * rho * cd * chord * relativeSpeed * relativeSpeed;


        //Open foam code is not decomposing Lift and drag
        //from wind axis to blade axis
        //Project over normal components
        vector normalForce(0,0,lift * cos(phi) - drag * sin(phi));
        vector tangentialForce(lift*sin(phi) + drag * cos(phi),0,0);
        
        pressOnPoints[i] = normalForce + tangentialForce;
        //Back to global ref frame
        pressOnPoints[i] = transform(bladeTensor,pressOnPoints[i]);

        //result.force += totalAerForce;
        //result.torque += (tangentialForce.x() * radius);

        //Back to global ref frame
        //totalAerForce = transform(bladeTensor,totalAerForce);

        //Add source term
        //force[celli] = -totalAerForce / volume;
        
    }

    forAll(rotorMesh_->rotorCells(),i)
    {
        const auto & rCell = rotorMesh_->rotorCells()[i];
        const auto & triangles = rCell.tri();
        vector bladeForce{0,0,0};
        label celli = rotorMesh_->cells()[i];
        forAll(triangles,j)
        {
            bladeForce += rCell.triArea()[j] * 
            (pressOnPoints[triangles[j][0]] + 
             pressOnPoints[triangles[j][1]] + 
             pressOnPoints[triangles[j][2]]) / 3;
        }
        bladeForce = pressOnPoints[rCell.center()] * rCell.area();
        result.force += bladeForce;

        result.torque += bladeForce.x() * cylPoints[rCell.center()].x();
        
        force[celli] = -bladeForce/rCell.volume();


        radiusField[celli] = cylPoints[rCell.center()].x();
    }
    Info<<"FORCE: "<<result.force<<endl;

    if(rotorMesh_->mesh().time().writeTime())
    {
        aoaField.write();
        clField.write();
        cdField.write();
        radiusField.write();
    }

    result.power = result.torque * omega;

    result.updateEta(this->refV);
    //Use J definition Vref/(rps * D)
    result.updateJ(this->refV,omega,rotorDiscrete_.geometry().radius());
    result.updateCT(this->refRho,omega,rotorDiscrete_.geometry().radius());
    result.updateCP(this->refRho,omega,rotorDiscrete_.geometry().radius());

    return result;

}
