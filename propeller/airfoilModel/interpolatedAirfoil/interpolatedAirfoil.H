#ifndef PROPELLER_AIRFOIL_MODEL_INTERPOLATED_AIRFOIL
#define PROPELLER_AIRFOIL_MODEL_INTERPOLATED_AIRFOIL


#include "airfoilModel.H"
#include "airfoilModelList.H"
#include "scalar.H"
#include "dictionary.H"

namespace Foam
{


class interpolatedAirfoil
    :  public interpolated<scalar,const airfoilModel*>, public airfoilModel
{
public:
    TypeName("interpolatedAirfoil");

    interpolatedAirfoil();
    interpolatedAirfoil(word name, const dictionary& dict);
    interpolatedAirfoil(interpolated<scalar,const airfoilModel*>& other);

    void build(const airfoilModelList& airfoilList);

    scalar cl(scalar alfaRad, scalar reynolds, scalar mach) const; 
    scalar cd(scalar alfaRad, scalar reynolds, scalar mach) const;

private:
    List<word> airfoilNames_;
};


}


#endif /* PROPELLER_AIRFOIL_MODEL_INTERPOLATED_AIRFOIL */

