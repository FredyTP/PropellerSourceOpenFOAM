#include "fileSampler.H"

#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "csvTable.H"
#include "ClosestNeighbor.H"

namespace Foam
{
 
defineTemplateTypeNameWithName(fileSampler<scalar>,"fileSampler");
addTemplatedToRunTimeSelectionTable(diskSampler,fileSampler,scalar,dictionary);

defineTemplateTypeNameWithName(fileSampler<vector>,"fileSampler");
addTemplatedToRunTimeSelectionTable(diskSampler,fileSampler,vector,dictionary);


template<class fType>
fileSampler<fType>::fileSampler(const dictionary& dict,const rotorGrid* rGrid,const rotorFvMeshSel* rMesh)
    : domainSampler<fType>(dict,rGrid,rMesh)
{
    this->read(dict);
}

template<>
bool fileSampler<vector>::read(const dictionary &dict)
{
    Info.stream().incrIndent();


    fileName csvpath = dict.getOrDefault<fileName>("csv","");
    if(csvpath!="")
    {
        csvTable<scalar,word> table(true);
        table.readFile(csvpath);

        List<scalar> vx = table.col("VelocityX");
        List<scalar> vy = table.col("VelocityY");
        List<scalar> vz = table.col("VelocityZ");

        List<scalar> x = table.col("X");
        List<scalar> y = table.col("Y");
        List<scalar> z = table.col("Z");


        List<vector> vel(vx.size());
        List<FixedList<scalar,3>> pos(vx.size());
        forAll(vel,i)
        {
            vel[i]=vector({vx[i],vy[i],vz[i]});
            pos[i]=FixedList<scalar,3>({x[i],y[i],z[i]});

        }


        ClosestNeighbor<scalar,vector,3> interp(pos,vel);

        const List<vector>& cellCenter = this->rMesh_->mesh().C();
        volVectorField U
        (
            IOobject
            (
                "Uread",
                this->rMesh_->mesh().time().timeName(),
                this->rMesh_->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->rMesh_->mesh(),
            dimensionedVector(dimVelocity,Zero)
        );
            


        forAll(cellCenter,j)
        {
            U[j] =interp.interpolate(FixedList<scalar,3>({cellCenter[j].x(),cellCenter[j].y(),cellCenter[j].z()})).value();
        }

        U.write();

        uFromFile_ = autoPtr<volVectorField>::New(std::move(U));

    }
    else
    {
        fileName file = dict.get<fileName>("file");

        volVectorField U
        (
            IOobject
            (
                file,
                this->rMesh_->mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            this->rMesh_->mesh()
        );
        uFromFile_ = autoPtr<volVectorField>::New(std::move(U));
    }
    
    return true;
    
}

template<>
bool fileSampler<scalar>::read(const dictionary &dict)
{
    Info.stream().incrIndent();


    fileName csvpath = dict.getOrDefault<fileName>("csv","");
    if(csvpath!="")
    {
        csvTable<scalar,word> table(true);
        table.readFile(csvpath);

        List<scalar> vx = table.col("scalar");

        List<scalar> x = table.col("X");
        List<scalar> y = table.col("Y");
        List<scalar> z = table.col("Z");


        List<FixedList<scalar,3>> pos(vx.size());
        forAll(pos,i)
        {
            pos[i]=FixedList<scalar,3>({x[i],y[i],z[i]});

        }


        ClosestNeighbor<scalar,scalar,3> interp(pos,vx);

        const List<vector>& cellCenter = this->rMesh_->mesh().C();
        volScalarField U
        (
            IOobject
            (
                "scalarRead",
                this->rMesh_->mesh().time().timeName(),
                this->rMesh_->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->rMesh_->mesh(),
            dimensionedScalar(dimless,Zero)
        );

        forAll(cellCenter,j)
        {
            U[j] = interp.interpolate(FixedList<scalar,3>({cellCenter[j].x(),cellCenter[j].y(),cellCenter[j].z()})).value();
        }

        U.write();

        uFromFile_ = autoPtr<volScalarField>::New(std::move(U));

    }
    else
    {
        fileName file = dict.get<fileName>("file");

        volScalarField U
        (
            IOobject
            (
                file,
                this->rMesh_->mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            this->rMesh_->mesh()
        );
        uFromFile_ = autoPtr<volScalarField>::New(std::move(U));
    }
    
    return true;
    
}

template<class fType>
const Field<fType>& fileSampler<fType>::sampleField(const GeometricField<fType, fvPatchField, volMesh>& U)
{
    return domainSampler<fType>::sampleField(*uFromFile_);
}

}

