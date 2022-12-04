

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
    propellerModel_(),
    rotorMesh_(mesh_)
{
    this->read(dict);
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
        

        //- Read geometry data
        rotorGeometry_.radius = coeffs_.getOrDefault<scalar>("radius", NO_RADIUS);
        rotorGeometry_.direction = coeffs_.getOrDefault("direction",vector(Zero));
        rotorGeometry_.center = coeffs_.getOrDefault("center",vector(Zero));

        //- Read propeller Model
        const dictionary& propellerModelDict = coeffs_.subDict("propellerModel");
        propellerModel_ = propellerModel::New(propellerModelDict);

        //- If no radius from dict, check on rotor Model
        if(rotorGeometry_.radius == NO_RADIUS)
        {
            rotorGeometry_.radius = propellerModel_->radius();
        }

        //READ ROTOR MESH
        const dictionary& rotorMeshDict = coeffs_.subDict("rotorMesh");
        rotorMesh_.read(rotorMeshDict);
        rotorMesh_.build(rotorGeometry_);
        //- Building rotor mesh may or may not modify rotorGeometry
        
        propellerModel_->setRotorMesh(&rotorMesh_);
        propellerModel_->build(rotorGeometry_);

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

    propellerModel_->calculate(force);

    eqn+=force;

    //If its time to write into files
    if(mesh_.time().writeTime())
    {
        //To save force field into a dict
        force.write();
    }
    
}