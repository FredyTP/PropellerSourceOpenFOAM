#ifndef PROPELLER_AIRFOIL_MODEL_CYLINDER_SIMPLE_AIRFOIL
#define PROPELLER_AIRFOIL_MODEL_CYLINDER_SIMPLE_AIRFOIL
#ifndef PROPELLER_AIRFOIL_MODEL_SIMPLE_AIRFOIL_SIMPLE_AIRFOIL
#define PROPELLER_AIRFOIL_MODEL_SIMPLE_AIRFOIL_SIMPLE_AIRFOIL

#include "airfoilModel.H"
#include "scalar.H"
#include "dictionary.H"

namespace Foam
{




class cylinderAirfoil
: public airfoilModel
{
public:
    
    TypeName("cylinderAirfoil");


    cylinderAirfoil(word name, scalar cd0);
    cylinderAirfoil(word name, const dictionary& dict);

    virtual scalar cl(scalar alfaRad, scalar reynolds, scalar mach) const;
    virtual scalar cd(scalar alfaRad, scalar reynolds, scalar mach) const;

    bool read(const dictionary& dict);
private:
    scalar cd0_;

};


}


#endif /* PROPELLER_AIRFOIL_MODEL_SIMPLE_AIRFOIL_SIMPLE_AIRFOIL */


#endif /* PROPELLER_AIRFOIL_MODEL_CYLINDER_SIMPLE_AIRFOIL */
