#ifndef PROPELLER_AIRFOIL_MODEL_INTERPOLATED_AIRFOIL_INTERPOLATED_AIRFOIL
#define PROPELLER_AIRFOIL_MODEL_INTERPOLATED_AIRFOIL_INTERPOLATED_AIRFOIL


#include "airfoilModel.H"
#include "airfoilModelList.H"
#include "scalar.H"
#include "dictionary.H"

namespace Foam
{

/**
 * @brief Airfoil class made out of interpolated airfoils
 * 
 */
class interpolatedAirfoil
    :  public Interpolated<scalar,const airfoilModel*>, public airfoilModel
{
public:
    TypeName("interpolatedAirfoil");

    interpolatedAirfoil();
    interpolatedAirfoil(word name, const dictionary& dict);
    interpolatedAirfoil(Interpolated<scalar,const airfoilModel*>& other);

    void build(const airfoilModelList& airfoilList);

    scalar cl(scalar alfaRad, scalar reynolds, scalar mach) const; 
    scalar cd(scalar alfaRad, scalar reynolds, scalar mach) const;

private:
    List<word> airfoilNames_;
};


}


#endif /* PROPELLER_AIRFOIL_MODEL_INTERPOLATED_AIRFOIL_INTERPOLATED_AIRFOIL */

