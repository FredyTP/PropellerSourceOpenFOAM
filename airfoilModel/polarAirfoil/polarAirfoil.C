#include "polarAirfoil.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"

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

bool polarAirfoil::read(const dictionary& dict)
{

    Info<<"Reading polar airfoil data for: " << this->airfoilName() << endl;

    bool ok=true;
    ok &= dict.readEntry("file",file_);

    // Filename  - ( Re - Ma)
    List<Tuple2<word,FixedList<scalar,2>>> polarFiles;

    //Read airfoil data
    if(!file_.empty())
    {
        Info<<"Reading polar data from: "<<file_<<endl;
        Foam::IFstream is(file_);
        is  >> polarFiles;
    }

    //Build airfoil polars
    forAll(polarFiles,i)
    {
        auto polarf = polarFiles[i];
        fileName polarfile = polarf.first();
        scalar Re = polarf.second()[0];
        scalar Ma = polarf.second()[1];

        fileName polarpath = file_;
        polarpath.replace_name(polarfile);
        
        polars_.push_back(new polar("lineal",polarpath,Re,Ma));
    }

    return ok;
}
scalar polarAirfoil::cl(scalar alfaRad, scalar reynolds, scalar mach) const
{
    if(polars_.size>()0)
    {
        return polars_[0]->cl(alfaRad);
    }
    return 0;

}
scalar polarAirfoil::cd(scalar alfaRad, scalar reynolds, scalar mach) const 
{
    if(polars_.size()>0)
    {
        return polars_[0]->cd(alfaRad);
    }
    return 0;
}


}