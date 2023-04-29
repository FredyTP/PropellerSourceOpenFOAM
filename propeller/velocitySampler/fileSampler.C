#include "fileSampler.H"

#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "csvTable.H"
#include "closestNeighbor.H"

namespace Foam
{
 
    defineTypeNameAndDebug(fileSampler,0);
    addToRunTimeSelectionTable(velocitySampler,fileSampler, dictionary);


fileSampler::fileSampler(const dictionary& dict,const rotorDiscrete* rDiscrete_,const rotorFvMeshSel* rMesh_)
    : offsetSampler(dict,rDiscrete_,rMesh_)
{
    this->read(dict);
}

bool fileSampler::read(const dictionary &dict)
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


        closestNeighbor<scalar,vector,3> interp(pos,vel);

        const List<vector>& cellCenter = velocitySampler::rMesh->mesh().C();
        volVectorField U
        (
            IOobject
            (
                "Uread",
                velocitySampler::rMesh->mesh().time().timeName(),
                velocitySampler::rMesh->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            velocitySampler::rMesh->mesh(),
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
                velocitySampler::rMesh->mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            velocitySampler::rMesh->mesh()
        );
        uFromFile_ = autoPtr<volVectorField>::New(std::move(U));
    }
    
    return true;
    
}

const vectorField& fileSampler::sampleVelocity(const volVectorField& U)
{
    return offsetSampler::sampleVelocity(*uFromFile_);
}

}
