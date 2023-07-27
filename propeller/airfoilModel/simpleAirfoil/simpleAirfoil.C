#include "simpleAirfoil.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    //set and define the type name "typeName"
    defineTypeNameAndDebug(simpleAirfoil,0);

    //Add to run time table to dynamically select the class based on
    //the dictionary airfoilModel atribute

    addToRunTimeSelectionTable(airfoilModel,simpleAirfoil,dictionary);

simpleAirfoil::simpleAirfoil(word name, scalar cl0,scalar cl_max,scalar cd_max,scalar cd0)
:   airfoilModel(name),
    cl0_(cl0),
    cl_max_(cl_max),
    cd_max_(cd_max),
    cd0_(cd0)
{

}
simpleAirfoil::simpleAirfoil(word name, const dictionary& dict)
:   airfoilModel(name,dict)
{
    this->read(dict);
}

bool simpleAirfoil::read(const dictionary& dict)
{
    Info<<"Reading simple airfoil data for:" << this->airfoilName() << endl;

    bool ok=true;
    ok &= dict.readEntry("cl0",cl0_);
    ok &= dict.readEntry("cd0",cd0_);
    ok &= dict.readEntry("cl_max",cl_max_);
    ok &= dict.readEntry("cd_max",cd_max_);

    return ok;
}
scalar simpleAirfoil::cl(scalar alfaRad, scalar reynolds, scalar mach) const
{
    return cl0_ + cl_max_ * sin(alfaRad)*cos(alfaRad);
}
scalar simpleAirfoil::cd(scalar alfaRad, scalar reynolds, scalar mach) const 
{
    return cd0_ + cd_max_ * pow(sin(alfaRad),2);
}


}


