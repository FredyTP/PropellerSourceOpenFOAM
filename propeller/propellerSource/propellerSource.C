

#include "propellerSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"
#include "propellerResult.H"

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
    rotorFvMeshSel_(mesh),
    force(
        IOobject
        (
            fv::option::name_ + ":rotorForce",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector(dimVelocity/dimTime, Zero)
    ),
    propResult_(name,mesh)
{
    //Output propeller adimensional parameter definition
    propellerResult::OutputDefinition(Info)<<endl;
    this->read(dict);
}


bool Foam::fv::propellerSource::read(const dictionary& dict)
{
    //Reads fv::option and saves propeller content in Coeffs_ dict
    if(fv::option::read(dict))
    {
        Info<<"Reading propeller source"<<endl;

        /*---------FV OPTIONS RELATED STUFF--------*/
        fv::option::resetApplied();

        //Read fields names to apply the source, if not present
        //source won't be apllied
        coeffs_.readEntry("fields", fieldNames_);
        Info<<"Applying to fields: "<<fieldNames_<<endl;
        fv::option::resetApplied();
        /*-----------------------------------------*/

        kDot_ = coeffs_.getOrDefault("kDot",0);

        /*----------READ USER SPECIFIED ROTOR GEOMETRY---------------*/
        //- Read geometry data if present
        rotorGeometry_.readIfPresent(coeffs_.subOrEmptyDict("geometry"));

        /*----------READ USER DESIRED ROTOR MODEL(BEMT ...)---------------*/
        //- Read propeller Model
        propellerModel_ = propellerModel::New(coeffs_.subDict("propellerModel"));

        /*----------READ FV ROTOR MESH CONFIG---------------*/
        rotorFvMeshSel_.read(coeffs_.subDict("rotorMesh"));

        /*-----BUILD MESH AND UPDATE PROVIDED ROTOR GEOMETRY-----*/
        rotorFvMeshSel_.build(rotorGeometry_);  //- Building rotor mesh may or may not modify rotorGeometry
       
        /*-----CONFIGURE THE PROPELLER MODEL-----*/
        propellerModel_->setRotorMesh(&rotorFvMeshSel_);
        //Build propeller model with 100% definitive rotorGeometry
        propellerModel_->build(rotorGeometry_);
        //Set reference properties for aerodynamic forces and adim variables
        propellerModel_->setRefRho(coeffs_.getOrDefault<scalar>("refRho",1.0));
        propellerModel_->setRefV(coeffs_.getOrDefault<scalar>("refV",1.0));
        //propellerModel_->rDiscrete().writeArea(this->name(),mesh_);
        Info<<endl;
        if(!rotorGeometry_.isReady())
        {
            FatalErrorInFunction
                <<"Rotor geometry data of "<< fv::option::name() <<" is not determined:"<<endl
                << rotorGeometry_
                << exit(FatalError);
        }

        Info<<"Rotor Geometry of: "<< fv::option::name()<<endl
        << rotorGeometry_;

        /*-----CREATE VELOCITY SAMPLING METHOD-----*/
        //Create velocitySampler for specified rotor discrete and mesh
        velSampler_ = velocitySampler::New(
                   dict.subDict("velocitySampler"),
                   &propellerModel_->rDiscrete(), 
                   &rotorFvMeshSel_);
        velSampler_->writeSampled(fv::option::name_); //Write to file sampling location

        
        /*----CREATE ROTOR DYNAMICS----*/
        dynamics_ = rotorDynamics::New(dict.subDict("dynamics"));

        return true;
    }

    return false;

}

void Foam::fv::propellerSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    Info<<"ScalarFieldi: "<<fieldi<<endl;
    Info<< name() << ": applying source to " << eqn.psi().name() << endl;
    /*volScalarField k
        (
            IOobject
            (
                fv::option::name_ + ":rotorForce",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimArea/pow(dimTime,3), Zero)
        );
    forAll(propellerModel_->rDiscrete().rotorCells(),i)
    {
        k[propellerModel_->rDiscrete().rotorCells()[i].celli()]=kDot_;
    }
    eqn+=k;*/

}
void Foam::fv::propellerSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    //Create force vector field from mesh and set dimensions

    Info<<"VectorFieldi: "<<fieldi<<endl;
    Info<< name() << ": applying source to " << eqn.psi().name() << endl;

    const volVectorField& Uin(eqn.psi());

    propResult_ = propellerModel_->calculate
        (
            velSampler_->sampleVelocity(Uin),
            dynamics_->angularVelocity(),
            force
        );

    dynamics_->integrate(mag(propResult_.torque),mesh_.time().deltaTValue());
    
    Info<<fv::option::name_<<": step parameters"<<endl;
    Info<<propResult_<<endl;
    //Add source term to the equation
    eqn+=force;

    //If its time to write into files
    if(mesh_.time().writeTime())
    {
        //To save force field into a dict
        //force.write();
        //Check if want to save sampled cells
        velSampler_->writeSampled(this->name());

    }

    
}
