

#include "propellerSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"

#include "cellSet.H"
#include "closestNeighbor.H"
// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(propellerSource,0);
        addToRunTimeSelectionTable(option,propellerSource,dictionary);
    }
}

const Foam::Enum
<
    Foam::fv::propellerSource::propellerModelType
>
Foam::fv::propellerSource::propellerModelTypeNames_
({
    {propellerModelType::pmActuatorDisk, "actuatorDisk"},
    {propellerModelType::pmBladeElement, "bladeElement"},
    {propellerModelType::pmVortexLaticce, "vortexLaticce"},
});

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::propellerSource::propellerSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(name,modelType,dict,mesh),
    propModel_(propellerModel::New(coeffs_)),
    rotorGeom_(autoPtr<rotorGeometry>::New(mesh,coeffs_)),
    airfoils_(coeffs_.subDict("airfoils"))
{
    read(dict);
}


bool Foam::fv::propellerSource::read(const dictionary& dict)
{
    //Reads fv::option and saves propeller content in Coeffs_ dict
    if(fv::option::read(dict))
    {
        std::cout<<"Reading propeller source"<<std::endl;
        fv::option::resetApplied();

        //Read fields names to apply the source, if not present
        //source won't be apllied
        coeffs_.readEntry("fields", fieldNames_);
        fv::option::resetApplied();
        //propModel_->read(Coeffs_);


        

        return true;
    }

    return false;

}

void Foam::fv::propellerSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    //Create force vector field from mesh and set dimensions
    volVectorField force
    (
        IOobject
        (
            name_ + ":rotorForce",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(eqn.dimensions()/dimVolume, Zero)
    );

    const scalarField& cellVolume = mesh_.V();
    const labelList& rotorCells = rotorGeom_->cells();
    forAll(rotorCells,i)
    {
        label celli = rotorCells[i];
        //force[celli]=-10000*rotorGeom_->direction();
        force[celli]=vector(0,-1000,0);
    }
    eqn-=force;

    //If its time to write into files
    if(mesh_.time().writeTime())
    {
        //To save force field into a dict
        force.write();
    }
    
}