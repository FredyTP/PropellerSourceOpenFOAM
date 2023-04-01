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
    
    integrationOrder = dict.getOrDefault<word>("integration","tri");
    correctCenters = dict.getOrDefault<bool>("correctCenters",false);
    refinementLevel = dict.getOrDefault<label>("borderRefinement",0);
}

Foam::scalar Foam::bladeElementModel::radius() const
{
    return bladeModel_.maxRadius();
}

void Foam::bladeElementModel::build(const rotorGeometry& rotorGeometry)
{
    rotorDiscrete_.buildCoordinateSystem(rotorGeometry);
    rotorDiscrete_.fromRotorMesh(*rotorFvMeshSel_,integrationOrder,correctCenters,refinementLevel);

    bladeModel_.setMaxRadius(rotorGeometry.radius());
}


Foam::propellerResult Foam::bladeElementModel::calculate(const vectorField& U,scalar angularVelocity, volVectorField& force)
{
    propellerResult result;
    //Puntos de la discretizacion
    const List<vector>& cylPoints = rotorDiscrete_.cylPoints();
    //Tensor de cada punto local to global
    const List<tensor>& bladeCS = rotorDiscrete_.localBladeCS();

    //Velocidad angular
    const scalar pi = Foam::constant::mathematical::pi;
    const scalar omega = angularVelocity;
    
    List<scalar> aoaList(cylPoints.size());
    List<vector> pressOnPoints(cylPoints.size(),vector(0,0,0));

    volScalarField radiusField
    (
        IOobject
        (
            "propeller:radius",
            rotorFvMeshSel_->mesh().time().timeName(),
            rotorFvMeshSel_->mesh()
        ),
        rotorFvMeshSel_->mesh(),
        dimensionedScalar(dimless, Zero)
    );
    volScalarField areaField
    (
        IOobject
        (
            "propeller:area",
            rotorFvMeshSel_->mesh().time().timeName(),
            rotorFvMeshSel_->mesh()
        ),
        rotorFvMeshSel_->mesh(),
        dimensionedScalar(dimArea, Zero)
    );

    forAll(cylPoints, i)
    {
        if(!rotorDiscrete_.integrationPoints()[i])
        {
            //If point is not need for integration, continue
            continue;
        }
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
    }



    const PtrList<rotorCell>& rotorCells = rotorDiscrete_.rotorCells();
    const scalarField& cellVol = rotorFvMeshSel_->mesh().V();

    forAll(rotorCells,i)
    {
        const auto & rCell = rotorCells[i];
        vector bladeForce = rCell.integrateField(pressOnPoints);
        label celli = rCell.celli();

        result.force += bladeForce;

        result.torque += bladeForce ^(cylPoints[rCell.center()].x() * bladeCS[rCell.center()].col<1>()) ;
        
        force[celli] = - bladeForce/cellVol[celli];


        radiusField[celli] = cylPoints[rCell.center()].x();
        areaField[celli] = rCell.area();
    }
    result.force = rotorDiscrete_.cartesian().localVector(result.force);
    result.torque = rotorDiscrete_.cartesian().localVector(result.torque);

    if(rotorFvMeshSel_->mesh().time().writeTime())
    {
        radiusField.write();
        areaField.write();
    }

    result.power = result.torque.z() * omega;

    result.updateEta(this->refV);
    //Use J definition Vref/(rps * D)
    result.updateJ(this->refV,omega,rotorDiscrete_.geometry().radius());
    result.updateCT(this->refRho,omega,rotorDiscrete_.geometry().radius());
    result.updateCP(this->refRho,omega,rotorDiscrete_.geometry().radius());

    return result;

}
