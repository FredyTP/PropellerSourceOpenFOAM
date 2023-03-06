#include "rotorDiscrete.H"
#include "cellSet.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "unitConversion.H"
#include "rotorGeometry.H"
namespace Foam
{

defineTypeNameAndDebug(rotorDiscrete,0);



rotorDiscrete::rotorDiscrete()
{

}
void rotorDiscrete::buildCoordinateSystem(const rotorGeometry& geometry)
{
    rotorGeometry_ = geometry;

    cylCS_ = coordSystem::cylindrical
                (
                    rotorGeometry_.center, //centerd to local
                    rotorGeometry_.direction, //z-axis 
                    rotorGeometry_.psiRef  //x-axis
                );
}
tensor rotorDiscrete::bladeLocalFromPoint(const point &localPoint) const
{
    //z- up, y -outwards from center, x perpendicular y,z (leading edge to trailing edge)
    point global,origin;

    global = cylCS_.globalPosition(localPoint);
    origin = cylCS_.origin();
    
    tensor rotTensor(cylCS_.R());

    //z-axis
    const vector ax3 = rotTensor.col<2>(); // == e3 (already normalized)

    //y-axis (radial direction)
    vector ax2(global - origin);

    ax2.removeCollinear(ax3);

    const scalar magAxis2(mag(ax2));

    // Trap zero size and colinearity
    if (magAxis2 < SMALL)
    {
        return rotTensor;
    }

    ax2 /= magAxis2;  // normalise

    // Replace with updated local axes

    rotTensor.col<0>(ax2^ax3);
    rotTensor.col<1>(ax2);

    return rotTensor;
}
void rotorDiscrete::fromRotorMesh(const rotorMesh &rotorMesh)
{
    Info<< "Building rotor Discrete from mesh" <<endl;
    discreteMode_ = discreteMode::dmMesh;
    const labelList cells = rotorMesh.cells();
    const auto &mesh = rotorMesh.mesh();
    
    cylPoints_.resize(cells.size());
    localBlade_.resize(cells.size());


    volScalarField volume
    (
        IOobject
        (
            "propeller:rotorVolume",
            rotorMesh.mesh().time().timeName(),
            rotorMesh.mesh()
        ),
        rotorMesh.mesh(),
        dimensionedScalar(dimVolume, Zero)
    );

    forAll(cells, i)
    {
        label celli = cells[i];
        vector cellCenter = mesh.C()[celli];

        //Project points over the disk
        //From global to local cyl position
        // Global -> local cyl
        cylPoints_[i] = cylCS_.localPosition(cellCenter);
        cylPoints_[i].z()=0;
        
        //Not the most efficient way to compute global coords
        localBlade_[i] = this->bladeLocalFromPoint(cylPoints_[i]);
        volume[celli] = mesh.V()[celli];
    }

    volume.write();
}
bool rotorDiscrete::read(const dictionary &dict)
{
    return false;
}

}