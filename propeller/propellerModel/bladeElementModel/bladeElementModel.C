#include "bladeElementModel.H"
#include "bladeSection.H"
#include "addToRunTimeSelectionTable.H"
#include "rotorGrid.H"

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
    bladeModel_(airfoils_,dict.subDict("bladeTest")),
    rotorDiscrete_(dict.subOrEmptyDict("discrete"))
{    
    dict.readEntry("nBlades",nBlades_);
    tipFactor_ = dict.getOrDefault<scalar>("tipFactor",1);
}

Foam::scalar Foam::bladeElementModel::radius() const
{
    return bladeModel_.maxRadius();
}

void Foam::bladeElementModel::build(const rotorGeometry& rotorGeometry)
{
    rotorDiscrete_.buildCoordinateSystem(rotorGeometry);
    rotorDiscrete_.fromMesh(*rotorFvMeshSel_);
    bladeModel_.setMaxRadius(rotorGeometry.radius());

    bladeModel_.writeBlade(300,"blade.csv");
}

Foam::propellerResult Foam::bladeElementModel::calculate(const vectorField& U,scalar angularVelocity, volVectorField& force)
{
    propellerResult result;
    //Puntos de la discretizacion
    const List<vector>& cylPoints = rotorDiscrete_.cylPoints();
    //Tensor de cada punto local to global
    const List<tensor>& bladeCS = rotorDiscrete_.localBladeCS();
    const PtrList<gridCell>& cells = rotorDiscrete_.grid.cells();


    volScalarField radialGrid
    (
        IOobject
        (
            "radialGrid",
            rotorFvMeshSel_->mesh().time().timeName(),
            rotorFvMeshSel_->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rotorFvMeshSel_->mesh(),
        dimensionedScalar(dimless, -1)
    );
    volScalarField cdField
    (
        IOobject
        (
            "cd",
            rotorFvMeshSel_->mesh().time().timeName(),
            rotorFvMeshSel_->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rotorFvMeshSel_->mesh(),
        dimensionedScalar(dimless, 0)
    );
    volScalarField azimutalGrid
    (
        IOobject
        (
            "azimutalGrid",
            rotorFvMeshSel_->mesh().time().timeName(),
            rotorFvMeshSel_->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rotorFvMeshSel_->mesh(),
        dimensionedScalar(dimless, -1)
    );

    volScalarField aoaField
    (
        IOobject
        (
            "aoa",
            rotorFvMeshSel_->mesh().time().timeName(),
            rotorFvMeshSel_->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rotorFvMeshSel_->mesh(),
        dimensionedScalar(dimless, -1)
    );

    volScalarField clField
    (
        IOobject
        (
            "cl",
            rotorFvMeshSel_->mesh().time().timeName(),
            rotorFvMeshSel_->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rotorFvMeshSel_->mesh(),
        dimensionedScalar(dimless, 0)
    );

    volVectorField interpolatedVel
    (
        IOobject
        (
            "relativeVel",
            rotorFvMeshSel_->mesh().time().timeName(),
            rotorFvMeshSel_->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rotorFvMeshSel_->mesh(),
        dimensionedVector(dimless, Zero)
    );

    volScalarField cellvol
    (
        IOobject
        (
            "cellVolume",
            rotorFvMeshSel_->mesh().time().timeName(),
            rotorFvMeshSel_->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rotorFvMeshSel_->mesh(),
        dimensionedScalar(dimVolume, Zero)
    );

    volScalarField torqueM
    (
        IOobject
        (
            "torque",
            rotorFvMeshSel_->mesh().time().timeName(),
            rotorFvMeshSel_->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rotorFvMeshSel_->mesh(),
        dimensionedScalar(dimless, Zero)
    );

    volVectorField computedForce
    (
        IOobject
        (
            "calcForce",
            rotorFvMeshSel_->mesh().time().timeName(),
            rotorFvMeshSel_->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rotorFvMeshSel_->mesh(),
        dimensionedVector(dimless, Zero)
    );


    //Velocidad angular
    const scalar pi = Foam::constant::mathematical::pi;
    const scalar omega = angularVelocity;
    
    scalar aoaMax = -VGREAT;
    scalar aoaMin = VGREAT;

    const scalarField& cellVol = rotorFvMeshSel_->mesh().V();

    forAll(cellVol,i)
    {
        cellvol[i]=cellVol[i];
    }

    List<vector> pressOnPoints(cells.size(),vector(0,0,0));
    scalar forcetot = 0.0;
    scalar momenttot = 0.0;
    //---CALCULATE VALUE ON INTEGRATION POINTS---//
    forAll(cells, i)
    {
        //Get local radius
        scalar radius = cells[i].radius();
        scalar chord,twist,sweep;
        interpolatedAirfoil airfoil;
        bladeModel_.sectionAtRadius(radius,chord,twist,sweep,airfoil);

        if(chord == 0)
        {
            continue;
        }
        
        
        scalar average_fact = nBlades_ * cells[i].dt() / (constant::mathematical::twoPi);
                
        //Local rotation tensor
        tensor bladeTensor = rotorDiscrete_.bladeLocalFromPoint(cells[i].center());

        //Global coordinate vector
        vector airVel = U[i];
        vector localAirVel = invTransform(bladeTensor,airVel);
        
        //Get blade velocity
        vector relativeBladeVel(omega*radius,0,0);

        //Get relative air velocity
        vector relativeVel = localAirVel + relativeBladeVel;

        //y component is radial, thus "doesn't contribute to aerodinamic forces"
        //relativeVel.y()=0;
        scalar relativeSpeed = mag(relativeVel);

        //Airspeed angle (positive when speed is from "below" airfoil)
        scalar phi = atan2(-relativeVel.z(),sign(relativeVel.x())*sqrt(pow(relativeVel.x(),2)+pow(relativeVel.y(),2))); //as starccm+
        //scalar phi = atan2(-relativeVel.z(),relativeVel.x());
        //Info<<"Phi: "<<phi<<endl;

        //Angle of atack
        scalar AoA = twist - phi;

        if(AoA < aoaMin) aoaMin = AoA;
        if(AoA > aoaMax) aoaMax = AoA;


        scalar rho = this->refRho;
        scalar nu = 1e-5;
        scalar re = rho*relativeSpeed*chord/nu;
        scalar c = 345;
        scalar mach = relativeSpeed/c;

        scalar cl = airfoil.cl(AoA,re,mach);
        scalar cd = airfoil.cd(AoA,re,mach);

        //Add tip factor effect:
        if(radius/rotorDiscrete_.geometry().radius()>=tipFactor_)
        {
            cl=0.0;
        }
       
        //Calculate aerodinamic forces
        scalar lift = average_fact * 0.5 * rho * cl * chord * relativeSpeed * relativeSpeed * cells[i].dr();
        scalar drag = average_fact * 0.5 * rho * cd * chord * relativeSpeed * relativeSpeed * cells[i].dr();

        //Open foam code is not decomposing Lift and drag
        //from wind axis to blade axis
        //Project over normal components
        vector normalForce(0,0,lift * cos(phi) - drag * sin(phi));
        vector tangentialForce(lift*sin(phi) + drag * cos(phi),0,0);
        
        pressOnPoints[i] = normalForce + tangentialForce;
        //Back to global ref frame
        pressOnPoints[i] = transform(bladeTensor,pressOnPoints[i]);       
        result.force += normalForce;
        result.torque += vector(0,0,tangentialForce.x()*radius);    

        const List<label>& cellis = cells[i].cellis();
        const List<scalar>& we = cells[i].weights();

        forAll(cellis,k)
        {
            force[cellis[k]] = - we[k]*pressOnPoints[i]/cellVol[cellis[k]];

            /*if(mag(force[cellis[k]])>10000)
            {
                Info<<cellis.size()<<endl;
                Warning<<"Force too strong"<<endl;
                Info<<"Force: "<<mag(force[cellis[k]])<<endl;
                Info<<"Radius: "<<radius<<endl;
                Info<<"ZoneForce: "<<pressOnPoints[i]<<endl;
                Info<<"Cell vol: "<<cellVol[cellis[k]]<<endl;
                Info<<"C: "<<chord<<endl;
                Info<<"RelSpeed: "<<relativeSpeed<<endl;
                Info<<"CL: "<<cl<<endl;
                Info<<"CD: "<<cd<<endl;
                Info<<"Re: "<<re<<endl;
                Info<<"AoA: "<<AoA<<endl;
                Info<<"U: "<<airVel<<endl;
                Info<<"Dr: "<<cells[i].dr()<<endl;
                Info<<"Dt: "<<cells[i].dt()<<endl;

            }*/
            cdField[cellis[k]]=cd;
            torqueM[cellis[k]]=tangentialForce.x()*radius/cells[i].dr();
            azimutalGrid[cellis[k]]=cells[i].theta();
            radialGrid[cellis[k]]=cells[i].radius();
            aoaField[cellis[k]]=AoA;
            clField[cellis[k]]=cl;
            interpolatedVel[cellis[k]]=U[i];
            computedForce[cellis[k]]=pressOnPoints[i]/cells[i].dr();
        }

    }


    reduce(aoaMax,maxOp<scalar>());
    reduce(aoaMin,minOp<scalar>());

    reduce(result.force,sumOp<vector>());
    reduce(result.torque,sumOp<vector>());

    result.power = result.torque.z() * omega;

    result.updateEta(this->refV);
    //Use J definition Vref/(rps * D)
    result.updateJ(this->refV,omega,rotorDiscrete_.geometry().radius());
    result.updateCT(this->refRho,omega,rotorDiscrete_.geometry().radius());
    result.updateCP(this->refRho,omega,rotorDiscrete_.geometry().radius());
    Info<< "- Max AoA: "<<aoaMax * 180/pi <<"º"<<endl;
    Info<< "- Min AoA: "<<aoaMin * 180/pi <<"º"<<endl;

    if(rotorFvMeshSel_->mesh().time().writeTime())
    {
        clField.write();
        azimutalGrid.write();
        radialGrid.write();
        aoaField.write();
        cdField.write();
        interpolatedVel.write();
        computedForce.write();
        cellvol.write();
        torqueM.write();
    }
    return result;

}
/*
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
    
    scalar aoaMax = -VGREAT;
    scalar aoaMin = VGREAT;

    List<vector> pressOnPoints(cylPoints.size(),vector(0,0,0));

    //---CALCULATE VALUE ON INTEGRATION POINTS---//
    forAll(cylPoints, i)
    {
        if(!rotorDiscrete_.integrationPoints()[i])
        {
            //If point is not need for integration, continue
            continue;
        }
        //Get local radius
        scalar radius = cylPoints[i].x();
        scalar chord,twist,sweep;
        interpolatedAirfoil airfoil;
        bladeModel_.sectionAtRadius(radius,chord,twist,sweep,airfoil);

        if(chord == 0)
        {
            continue;
        }
        
        
        scalar average_fact = nBlades_ / (2 * pi * radius);
                
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

        if(AoA < aoaMin) aoaMin = AoA;
        if(AoA > aoaMax) aoaMax = AoA;


        scalar rho = this->refRho;
        scalar nu = 1e-5;
        scalar re = rho*relativeSpeed*chord/nu;
        scalar c = 345;
        scalar mach = relativeSpeed/c;

        scalar cl = airfoil.cl(AoA,re,mach);
        scalar cd = airfoil.cd(AoA,re,mach);

        //Add tip factor effect:
        if(radius/rotorDiscrete_.geometry().radius()>=tipFactor_)
        {
            cl=0.0;
        }
       
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

    reduce(aoaMax,maxOp<scalar>());
    reduce(aoaMin,minOp<scalar>());

    Info<< "- Max AoA: "<<aoaMax * 180/pi <<"º"<<endl;
    Info<< "- Min AoA: "<<aoaMin * 180/pi <<"º"<<endl;

    const PtrList<rotorCell>& rotorCells = rotorDiscrete_.rotorCells();
    const scalarField& cellVol = rotorFvMeshSel_->mesh().V();

    //-----INTEGRATE PRESSURE FIELD------//
    forAll(rotorCells,i)
    {
        const auto & rCell = rotorCells[i];
        vector bladeForce = rCell.integrateField(pressOnPoints);
        label celli = rCell.celli();
        force[celli] = - bladeForce/cellVol[celli];

        result.force += bladeForce;
        result.torque += bladeForce ^(cylPoints[rCell.center()].x() * bladeCS[rCell.center()].col<1>()) ;        
    }

    //To rotor local coordinates
    result.force = rotorDiscrete_.cartesian().localVector(result.force);
    result.torque = rotorDiscrete_.cartesian().localVector(result.torque);

    reduce(result.force,sumOp<vector>());
    reduce(result.torque,sumOp<vector>());

    result.power = result.torque.z() * omega;

    result.updateEta(this->refV);
    //Use J definition Vref/(rps * D)
    result.updateJ(this->refV,omega,rotorDiscrete_.geometry().radius());
    result.updateCT(this->refRho,omega,rotorDiscrete_.geometry().radius());
    result.updateCP(this->refRho,omega,rotorDiscrete_.geometry().radius());

    return result;

}*/
