#ifndef PROPELLER_AIRFOIL_MODEL_VKT_AIRFOIL_VKT_AIRFOIL
#define PROPELLER_AIRFOIL_MODEL_VKT_AIRFOIL_VKT_AIRFOIL

#include "airfoilModel.H"
#include "scalar.H"
#include "dictionary.H"

namespace Foam
{



/**
 * Semi-empirical analitic airfoil model over 360º range base on paper:
 * An analytical model for airfoil aerodynamic 
 * characteristics over the entire 360 degrees angle of attack range.
 * V.K. Truong
 * DOI: 10.1063/1.5126055
 * 
 * This model is valid for:
 * M<<1
 * Stationary (no account for unsteady effects like dynamic stall)
*/
class vktAirfoil
: public airfoilModel
{
public:
    
    TypeName("vktAirfoil");


    vktAirfoil(const word name, scalar cl0,scalar cl_alfa,scalar K,scalar cd0);
    vktAirfoil(const word name, const dictionary& dict);

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


#endif /* PROPELLER_AIRFOIL_MODEL_VKT_AIRFOIL_VKT_AIRFOIL */
