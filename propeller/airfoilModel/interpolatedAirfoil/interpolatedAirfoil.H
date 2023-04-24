#ifndef PROPELLER_AIRFOIL_MODEL_INTERPOLATED_AIRFOIL
#define PROPELLER_AIRFOIL_MODEL_INTERPOLATED_AIRFOIL


#include "airfoilModel.H"
#include "scalar.H"
#include "dictionary.H"

namespace Foam
{




class interpolatedAirfoil
    :  public interpolated<scalar,const airfoilModel*>
{
public:
    //TypeName("interpolatedAirfoil");

    interpolatedAirfoil();
    interpolatedAirfoil(interpolated<scalar,const airfoilModel*>& other);

    scalar cl(scalar alfaRad, scalar reynolds, scalar mach) const; 
    scalar cd(scalar alfaRad, scalar reynolds, scalar mach) const;
};


}


#endif /* PROPELLER_AIRFOIL_MODEL_INTERPOLATED_AIRFOIL */

