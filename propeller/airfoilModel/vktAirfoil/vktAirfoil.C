#include "vktAirfoil.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    //set and define the type name "typeName"
    defineTypeNameAndDebug(vktAirfoil,0);

    //Add to run time table to dynamically select the class based on
    //the dictionary airfoilModel atribute

    addToRunTimeSelectionTable(airfoilModel,vktAirfoil,dictionary);

vktAirfoil::vktAirfoil(const word name, scalar cl0, scalar cl_alfa, scalar K, scalar cd0)  
:   airfoilModel(name),
    cl0_(cl0),
    dcl_dalfa_(cl_alfa),
    K_(K),
    cd0_(cd0)
{

}
vktAirfoil::vktAirfoil(const word name, const dictionary& dict)
:   airfoilModel(name,dict)
{
    this->read(dict);
}

bool vktAirfoil::read(const dictionary& dict)
{

    Info<<"Reading simple airfoil data for:" << this->airfoilName() << endl;

    bool ok=true;
    ok &= dict.readEntry("cl0",cl0_);
    ok &= dict.readEntry("cd0",cd0_);
    ok &= dict.readEntry("cl_alfa",dcl_dalfa_);
    ok &= dict.readEntry("K",K_);

    return ok;
}
scalar vktAirfoil::cl(scalar alfaRad, scalar reynolds, scalar mach) const
{
    return cl0_ + dcl_dalfa_ * alfaRad;
}
scalar vktAirfoil::cd(scalar alfaRad, scalar reynolds, scalar mach) const
{
    scalar CL = cl(alfaRad,reynolds,mach);

    return cd0_ + K_ * CL* CL;
}


}