#include "polarAirfoil.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "linearInterpolation.H"

namespace Foam
{
    //set and define the type name "typeName"
    defineTypeNameAndDebug(polarAirfoil,0);

    //Add to run time table to dynamically select the class based on
    //the dictionary airfoilModel atribute

    addToRunTimeSelectionTable(airfoilModel,polarAirfoil,dictionary);

polarAirfoil::polarAirfoil(const word name, const dictionary& dict)
:   airfoilModel(name)
{
    this->read(dict);
}

/*
bool polarAirfoil::readTable(const dictionary& dict)
{
    Info<<"Reading polar airfoil data for: " << this->airfoilName() << endl;
    bool ok=true;
    ok &= dict.readEntry("file",file_);
    word extrapolation = dict.getOrDefault<word>("extrapolation","polar");
    //Read airfoil data
    if(!file_.empty())
    {
        Info<<"Reading polar data from: "<<file_<<endl;
        csvTable<scalar,word> csv(true);
        csv.readFile(file_);

        auto Re = csv.col2("Re");
        auto cl = csv.col("Cl");
        auto cd = csv.col("Cd");
        auto aoa = csv.col("AoA")
        
    }
}*/

bool polarAirfoil::read(const dictionary& dict)
{

    Info<<"Reading polar airfoil data for: " << this->airfoilName() << endl;

    bool ok=true;
    ok &= dict.readEntry("file",file_);
    word extrapolation = dict.getOrDefault<word>("extrapolation","polar");

    // Filename  - ( Re - Ma)
    List<Tuple2<word,FixedList<scalar,2>>> polarFiles;

    //Read airfoil data
    if(!file_.empty())
    {
        Info<<"Reading polar data from: "<<file_<<endl;
        Foam::IFstream is(file_);
        is  >> polarFiles;
    }
    polars_.resize(polarFiles.size());
    //Build airfoil polars
    forAll(polarFiles,i)
    {
        auto polarf = polarFiles[i];
        fileName polarfile = polarf.first();
        scalar Re = polarf.second()[0];
        scalar Ma = polarf.second()[1];

        fileName polarpath = file_;
        polarpath.replace_name(polarfile);
        auto ptrPolar = polar::New(extrapolation,"lineal",polarpath,Re,Ma);
        polars_[i].reset(ptrPolar.release());
    }

    //Create list to create interpolation table
    List<polar*> polarList(polars_.size());
    List<List<scalar>> ReMa(polars_.size());
    forAll(polars_,i)
    {
        polarList[i]=polars_[i].get();
        ReMa[i].setSize(2);
        ReMa[i][0]=polars_[i]->reynolds();
        ReMa[i][1]=polars_[i]->mach();
    }
    polarInterpolated = autoPtr<interpolationTable<scalar,polar*,2>>::NewFrom<linearInterpolation<scalar,polar*,2>>();
    polarInterpolated->setRawData(ReMa,polarList);

    return ok;
}
scalar polarAirfoil::cl(scalar alfaRad, scalar reynolds, scalar mach) const
{
    return polarInterpolated->interpolate({reynolds,mach}).value([=](scalar val,polar* p){return val*p->cl(alfaRad);});
}
scalar polarAirfoil::cd(scalar alfaRad, scalar reynolds, scalar mach) const
{
    return polarInterpolated->interpolate({reynolds,mach}).value([=](scalar val,polar* p){return val*p->cd(alfaRad);});
}


}