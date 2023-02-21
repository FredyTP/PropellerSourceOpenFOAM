#ifndef PROPELLER_AIRFOIL_MODEL_POLAR_AIRFOIL_POLAR_AIRFOIL
#define PROPELLER_AIRFOIL_MODEL_POLAR_AIRFOIL_POLAR_AIRFOIL

#include "airfoilModel.H"
#include "scalar.H"
#include "dictionary.H"

namespace Foam
{




class polarAirfoil
: public airfoilModel
{
public:
    
    TypeName("polarAirfoil");


    polarAirfoil(const word name, scalar cl0,scalar cl_alfa,scalar K,scalar cd0);
    polarAirfoil(const word name, const dictionary& dict);

    virtual scalar cl(scalar alfaRad, scalar reynolds, scalar mach) const ;
    virtual scalar cd(scalar alfaRad, scalar reynolds, scalar mach) const ;

    bool read(const dictionary& dict);
private:

    scalar cl0_;
    scalar dcl_dalfa_;
    scalar K_;
    scalar cd0_;

};


}


#endif /* PROPELLER_AIRFOIL_MODEL_POLAR_AIRFOIL_POLAR_AIRFOIL */
