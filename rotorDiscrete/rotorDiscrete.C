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

    coordSys_ = coordSystem::cylindrical
                (
                    rotorGeometry_.center, 
                    rotorGeometry_.direction,
                    rotorGeometry_.psiRef
                );
}
void rotorDiscrete::fromRotorMesh(const rotorMesh& rotorMesh)
{
    Info<< "Building rotor Discrete from mesh" <<endl;
    discreteMode_ = discreteMode::dmMesh;
    const labelList cells = rotorMesh.cells();
    const auto &mesh = rotorMesh.mesh();
    
    cylPoints_.resize(cells.size());


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
        cylPoints_[i] = coordSys_.localPosition(cellCenter);
        volume[celli] = mesh.V()[celli];
    }

    volume.write();


    


}
bool rotorDiscrete::read(const dictionary &dict)
{
    
}

}