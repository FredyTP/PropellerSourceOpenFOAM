

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
        fv::option::resetApplied();
        /*-----------------------------------------*/


        /*----------READ USER SPECIFIED ROTOR GEOMETRY---------------*/
        //- Read geometry data, SET TO EMPTY VALUE IF NOT PRESENT
        rotorGeometry_.readIfPresent(coeffs_);

        /*----------READ USER DESIRED ROTOR MODEL(BEMT ...)---------------*/
        //- Read propeller Model
        const dictionary& propellerModelDict = coeffs_.subDict("propellerModel");
        propellerModel_ = propellerModel::New(propellerModelDict);

        /*---------IF NO RADIUS SPECIFIED CHECK FROM MODEL DATA---------------*/
        //- If no radius from dict, check on rotor Model
        if(!rotorGeometry_.isRadiusSet())
        {
            //From blade data if not adimensional
            if(propellerModel_->radius() != NO_RADIUS) 
                rotorGeometry_.setRadius(propellerModel_->radius());
        }

        /*----------READ FV ROTOR MESH CONFIG---------------*/
        const dictionary& rotorMeshDict = coeffs_.subDict("rotorMesh");
        rotorMesh_.read(rotorMeshDict);

        /*-----BUILD MESH AND UPDATE PROVIDED ROTOR GEOMETRY-----*/
        rotorMesh_.build(rotorGeometry_);  //- Building rotor mesh may or may not modify rotorGeometry
       
        /*-----SET THE FV MESH TO THE ROTOR MODEL-----*/
        propellerModel_->setRotorMesh(&rotorMesh_);
        //Build propeller model with 100% definitive rotorGeometry
        propellerModel_->build(rotorGeometry_);

        Info<<"Rotor Geometry of: "<< this->name()<<endl
        << rotorGeometry_;

        if(!rotorGeometry_.isReady())
        {
            FatalErrorInFunction
                <<"Rotor geometry data of "<< this->name() <<" is not determined:"
                << exit(FatalError);
        }

        //Set reference properties for aerodynamic forces and adim variables
        propellerModel_->setRefRho(coeffs_.getOrDefault<scalar>("refRho",1.0));
        propellerModel_->setRefV(coeffs_.getOrDefault<scalar>("refV",1.0));

        /*-----CREATE VELOCITY SAMPLING METHOD-----*/
        //Create velocitySampler for specified rotor discrete and mesh
        velSampler_ = velocitySampler::New(
                   dict.subDict("velocitySampler"),
                   &propellerModel_->rDiscrete(), 
                   &rotorMesh_);
        velSampler_->writeSampled(name_); //Write to file sampling location

        
        /*----CREATE ROTOR DYNAMICS----*/
        dynamics_ = rotorDynamics::New(dict.subDict("dynamics"));



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

    const volVectorField& Uin(eqn.psi());

    propellerResult result;
    result = propellerModel_->calculate
        (
            velSampler_->sampleVelocity(Uin),
            dynamics_->angularVelocity(),
            force
        );

    dynamics_->integrate(result.torque,mesh_.time().deltaTValue());
    
    Info<<name_<<": step parameters"<<endl;
    Info<<result<<endl;
    //Add source term to the equation
    eqn+=force;

    //If its time to write into files
    if(mesh_.time().writeTime())
    {
        //To save force field into a dict
        force.write();
    }
    
}

