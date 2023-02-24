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


    localCS_ = coordSystem::cartesian
                (
                    rotorGeometry_.center,
                    rotorGeometry_.direction, //z-axis
                    rotorGeometry_.psiRef    //x-axis
                );

    cylCS_ = coordSystem::cylindrical
                (
                    vector(0,0,0), //centerd to local
                    vector(0,0,1), //z axis is the same
                    vector(1,0,0)  //x axis is the same
                );
}
void rotorDiscrete::fromRotorMesh(const rotorMesh& rotorMesh)
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
        //From global to local cyl position
        // Global -> local cart -> local cyl
        cylPoints_[i] = cylCS_.localPosition(localCS_.localPosition(cellCenter));
        cylPoints_[i].z()=0;
        
        localBlade_[i] = cylCS_.R(coordSystem::cylindrical::toCartesian(cylPoints_[i]));
        volume[celli] = mesh.V()[celli];
        Info<<cylPoints_[i]<<endl;
    }

    volume.write();
}
bool rotorDiscrete::read(const dictionary &dict)
{
    return false;
}

}