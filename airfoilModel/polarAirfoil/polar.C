#include "polar.H"
#include "linearInterpolation.H"
#include "csvTableReader.H"
#include "dictionary.H"
namespace Foam
{

polar::polar(const word interpolation, List<scalar> alpha, List<scalar> cl, List<scalar> cd, scalar Re, scalar Ma)
{
    cl_alpha =  autoPtr<regularInterpolation<1>>::NewFrom<linearInterpolation1>(alpha,cl);
    cd_alpha =  autoPtr<regularInterpolation<1>>::NewFrom<linearInterpolation1>(alpha,cd);
    reynolds_ = Re;
    mach_ = Ma;
}

polar::polar(const word interpolation,fileName filename, scalar Re, scalar Ma)
{

    //Config csv reader
    dictionary dict;
    dict.add("hasHeaderLine","true");
    dict.add("refColumn",0);
    dict.add("componentColumns","(0 1 2)");

    csvTableReader<vector> reader(dict);

    List<Tuple2<scalar,vector>> data;

    reader(filename,data);

    //File no found
    if(data.size()==0)
    {
        FatalErrorInFunction
                << "Cannot find file: " << filename
                << ". Required by polar class."<< nl
                << exit(FatalError);
    }

    //Data arrays
    List<scalar> alpha;
    List<scalar> cl;
    List<scalar> cd;

    alpha.resize(data.size());
    cl.resize(data.size());
    cd.resize(data.size());

    //Transfer data to arrays
    forAll(data,i)
    {
        vector vec = data[i].second();
        alpha[i]=vec[0];
        cl[i]=vec[1];
        cd[i]=vec[2];
    }

    //Create polar interpolations
    cl_alpha =  autoPtr<regularInterpolation<1>>::NewFrom<linearInterpolation1>(alpha,cl);
    cd_alpha =  autoPtr<regularInterpolation<1>>::NewFrom<linearInterpolation1>(alpha,cd);
}

scalar polar::cl(scalar alpha)
{
    return cl_alpha->interpolate(alpha);
}

scalar polar::cd(scalar alpha)
{
    return cd_alpha->interpolate(alpha);
}
polar::~polar()
{
    Info<<"Deleted polar"<<endl;
}

}
