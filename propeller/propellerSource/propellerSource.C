#include "propellerSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"
#include "propellerResult.H"
#include "cellSet.H"
#include "ClosestNeighbor.H"
#include "SquareMatrix.H"
#include "simpleMatrix.H"
#include "RegularInterpolation.H"
#include "csvTable.H"
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
        incrIndent(Info);
        /*---------FV OPTIONS RELATED STUFF--------*/
        fv::option::resetApplied();

        //Read fields names to apply the source, if not present
        //source won't be apllied
        coeffs_.readEntry("fields", fieldNames_);
        Info<<"Applying to fields: "<<fieldNames_<<endl;
        fv::option::resetApplied();
        /*-----------------------------------------*/

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

        /*-----BUILD PROPELLER MODEL AND COORDINATE SYSTEM-----*/
        rotorGeometry_.buildCS();
        propellerModel_->build(rotorGeometry_);

        /*-----CREATE VELOCITY SAMPLING METHOD-----*/
        //Create velocitySampler for specified rotor discrete and mesh
        Info<<endl;
        Info<<"Velocity sampler:"<<endl;
        velSampler_ = diskSampler<vector>::New(
                   dict.subDict("velocitySampler"),
                   propellerModel_->grid().get(), 
                   &rotorFvMeshSel_);

        velSampler_->writeSampled(fv::option::name_); //Write to file sampling location

        Info<<endl;
        Info<<"Density sampler (ignore if incompressible):"<<endl;
        dictionary densityDict;
        densityDict.add("type","domainSampler");
        if(dict.findDict("densitySampler"))
        {
            densityDict = dict.subDict("densitySampler");
        }
        densitySampler_ = diskSampler<scalar>::New(
                    densityDict,
                    propellerModel_->grid().get(), 
                    &rotorFvMeshSel_);
        
        decrIndent(Info);
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
    this->addSup(volScalarField::null(),eqn,fieldi);
}

void Foam::fv::propellerSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    force.ref(false)=dimensionedVector(dimVelocity/dimTime, Zero);

    Info<< name() << ": applying source to " << eqn.psi().name() << endl;
    Info<<fv::option::name_<<": step parameters"<<endl;

    const volVectorField& Uin(eqn.psi());

    if(propellerModel_->nextTimeStep(mesh_.time().deltaTValue()))
    {
        //If changes are made updating time step

        velSampler_->build();
        densitySampler_->build();
    }


    propResult_ = propellerModel_->calculate
        (
            velSampler_->sampleField(Uin),
            rho.empty() ? &densitySampler_->sampleField(rho) : nullptr,
            force
        );
    
    Info<<propResult_<<endl;
    //Add source term to the equation
    eqn+=force;

    //If its time to write into files
    if(mesh_.time().writeTime())
    {
        velSampler_->writeSampled(this->name());
    }

}
