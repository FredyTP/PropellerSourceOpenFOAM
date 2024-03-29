#include "cylinderAirfoil.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    //set and define the type name "typeName"
    defineTypeNameAndDebug(cylinderAirfoil,0);

    //Add to run time table to dynamically select the class based on
    //the dictionary airfoilModel atribute

    addToRunTimeSelectionTable(airfoilModel,cylinderAirfoil,dictionary);

cylinderAirfoil::cylinderAirfoil(word name, scalar cd0)  
:   airfoilModel(name),
    cd0_(cd0)
{

}
cylinderAirfoil::cylinderAirfoil(word name, const dictionary& dict)
:   airfoilModel(name,dict)
{
    this->read(dict);
}

bool cylinderAirfoil::read(const dictionary& dict)
{
    bool ok=true;
    ok &= dict.readEntry("cd0",cd0_);

    return ok;
}
scalar cylinderAirfoil::cl(scalar alfaRad, scalar reynolds, scalar mach) const
{
    return 0;
}
scalar cylinderAirfoil::cd(scalar alfaRad, scalar reynolds, scalar mach) const 
{
    return cd0_;
}


}

