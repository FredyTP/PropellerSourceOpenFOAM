#include "polar.H"
#include "linearInterpolation1D.H"

#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"
#include "Tuple2.H"
#include "vector.H"
#include "csvTableReader.H"
#include "mathematicalConstants.H"

namespace Foam
{
    defineTypeNameAndDebug(polar,0);

    //Define run time table for selecting derived types
    defineRunTimeSelectionTable(polar, dictionary);   

    addToRunTimeSelectionTable(polar,polar,dictionary);


polar::polar(const word interpolation, List<scalar>& alpha, List<scalar>& cl, List<scalar>& cd, scalar Re, scalar Ma)
{
    processData(alpha,cl,cd);

    cl_alpha =  autoPtr<regularInterpolation1D>::NewFrom<linearInterpolation1D>(alpha,cl);
    cd_alpha =  autoPtr<regularInterpolation1D>::NewFrom<linearInterpolation1D>(alpha,cd);
    reynolds_ = Re;
    mach_ = Ma;
}

polar::polar(const word interpolation,fileName filename, scalar Re, scalar Ma)
{
    reynolds_ = Re;
    mach_ = Ma;

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
        alpha[i]=vec[0]*constant::mathematical::pi/180;
        cl[i]=vec[1];
        cd[i]=vec[2];
    }

    processData(alpha,cl,cd);

    //Create polar interpolations
    cl_alpha =  autoPtr<regularInterpolation1D>::NewFrom<linearInterpolation1D>(alpha,cl);
    cd_alpha =  autoPtr<regularInterpolation1D>::NewFrom<linearInterpolation1D>(alpha,cd);
}

scalar polar::cl(scalar alpha)
{
    return cl_alpha->interpolate(alpha);
}

scalar polar::cd(scalar alpha)
{
    return cd_alpha->interpolate(alpha);
}
scalar polar::reynolds() const
{
    return reynolds_;
}
scalar polar::mach() const
{
    return mach_;
}
polar::~polar()
{
    Info<<"Deleted polar"<<endl;
}

void polar::processData(List<scalar> &alpha, List<scalar> &cl, List<scalar> &cd)
{
    //Check if alpha is in ascending order
    if(!std::is_sorted(alpha.begin(),alpha.end()))
    {
        //FatalError()
        return;
    }

    alpha_min = alpha[0];
    alpha_max = alpha[alpha.size()-1];

    cl_alpha_min=cl[0];
    cl_alpha_max=cl[cl.size()-1];
    cd_alpha_min=cd[0];
    cd_alpha_max=cd[cd.size()-1];

}

autoPtr<Foam::polar> polar::New(const word modelType, const word interpolation, const fileName filename, scalar Re, scalar Ma)
{
    Info<< "    Selecting " << typeName << " " << modelType << endl;

    //Find class contructor in tables
    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        /*FatalIOErrorInLookup
        (
            "polar",
            typeName,
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);*/
    }

    return autoPtr<Foam::polar>(ctorPtr(interpolation,filename,Re,Ma));
}

}
