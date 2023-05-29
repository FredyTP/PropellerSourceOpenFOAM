#include "propellerModel.H"

#include "axisAngleRotation.H"

namespace Foam
{

    defineTypeNameAndDebug(propellerModel,0);

    //Define run time table for selecting derived types
    defineRunTimeSelectionTable(propellerModel, dictionary);   

    bool propellerResult::definitionShown_  = false;


propellerModel::propellerModel
(
    const dictionary& dict,
    const word& name
)
{


}

tensor propellerModel::bladeTensor(const coordSystem::cylindrical &cylCS, const point &localPoint, scalar flapping, scalar sweep)
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
    return rotTensor.inner(rotX);
}

}