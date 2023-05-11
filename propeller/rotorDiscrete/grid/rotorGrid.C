#include "rotorGrid.H"
#include "polarGrid.H"
#include "bladeGrid.H"
#include "mathematicalConstants.H"
#include "regularInterpolation.H"

namespace Foam
{

rotorGrid::rotorGrid(const rotorGeometry& geometry, label nBlades)
{
    rotorGeometry_ = geometry;
    nBlades_ = nBlades;
    minRadius_ = geometry.innerRadius();
    maxRadius_ = geometry.radius();
}
void rotorGrid::setCenterFromClosestCell(const vectorField &cellCenter)
{
    forAll(cells_, i)
    {
        cells_[i].centerFromClosestCell(cellCenter, rotorGeometry_.cylindricalCS());
    }
}
tensor rotorGrid::bladeLocalFromPoint(const coordSystem::cartesian &carCS, const point &localPoint) 
{
    // z- up, y -outwards from center, x perpendicular y,z (leading edge to trailing edge)
    point global, origin;

    global = carCS.globalPosition(localPoint);
    origin = carCS.origin();

    tensor rotTensor(carCS.R());

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

    return rotTensor;
}
tensor rotorGrid::bladeLocalFromPoint(const coordSystem::cylindrical &cylCS, const point &localPoint) 
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

    return rotTensor;
}

autoPtr<rotorGrid> rotorGrid::New(const dictionary &dict, const rotorGeometry& geometry, scalar nBlades)
{
    word type = dict.get<word>("type");
    scalar innerRadius = geometry.innerRadius();
    scalar radius = geometry.radius();

    autoPtr<rotorGrid> ptr;
    if(type=="polarGrid")
    {
        Info<<"Creating polar grid"<<endl;
        scalar nRadial = dict.get<label>("nRadial");
        scalar nAzimutal = dict.get<label>("nAzimutal");
        ptr = autoPtr<rotorGrid>::NewFrom<polarGrid>(geometry,nBlades,nRadial,nAzimutal);
    }
    else if ( type == "meshGrid")
    {

    }
    else if (type == "bladeGrid")
    {
        Info<<"Creating bladeGrid"<<endl;
        scalar nRadial = dict.get<label>("nRadial");
        scalar nChord = dict.get<label>("nChord");
        scalar chord = 0.05;
        ptr = autoPtr<rotorGrid>::NewFrom<bladeGrid>(geometry,chord,nBlades,nRadial,nChord);
    }

    return ptr;
}
}

