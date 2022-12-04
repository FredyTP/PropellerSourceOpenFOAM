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
    const List<vector> cylPoints = rotorDiscrete_.cylPoints();
    double rpm = 10000;
    double omega = rpm*2*3.1416/60;
    forAll(cylPoints, i)
    {
        scalar radius = cylPoints[i].x();
        scalar bladeSpeed = omega*radius;
        scalar alpha = 10 *3.1416/180;
        scalar chord = bladeModel_.chordAtRadius(radius);
        scalar rho = 1.225;
        scalar nu = 1e-5;
        scalar re = rho*bladeSpeed*chord/nu;
        scalar c = 345;
        scalar mach = bladeSpeed/c;
        scalar cl = airfoils_.getAirfoil(0)->cl(alpha,re,mach);
        scalar cd = airfoils_.getAirfoil(0)->cd(alpha,re,mach);

        scalar area = rotorMesh_->areas()[i];
        label celli = rotorMesh_->cells()[i];
        scalar volume = rotorMesh_->mesh().V()[celli];
        scalar lift = 0.5 * rho * cl * chord * bladeSpeed * bladeSpeed * area /(2 * 3.1416 * radius);
        scalar drag = 0.5 * rho * cd * chord * bladeSpeed * bladeSpeed * area /(2 * 3.1416 * radius);
        
        totalLift += lift;
        totalDrag += drag;
    
        
        force[celli] = -lift * rotorDiscrete_.geometry().direction / volume;
    }

    Info<< "Total Lift: "<<totalLift<<endl;
    Info<< "Total Drag: "<<totalDrag<<endl;

}