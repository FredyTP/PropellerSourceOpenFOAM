#ifndef PROPELLER_AIRFOIL_MODEL_POLAR_AIRFOIL_POLAR_AIRFOIL
#define PROPELLER_AIRFOIL_MODEL_POLAR_AIRFOIL_POLAR_AIRFOIL

#include "airfoilModel.H"
#include "scalar.H"
#include "dictionary.H"
#include "polar.H"
#include "Tuple2.H"
#include "interpolationTable.H"
namespace Foam
{

class polarAirfoil
: public airfoilModel
{
public:
    
    TypeName("polarAirfoil");

    polarAirfoil(const word name, const dictionary& dict);

    virtual scalar cl(scalar alfaRad, scalar reynolds, scalar mach) const;
    virtual scalar cd(scalar alfaRad, scalar reynolds, scalar mach) const;

    bool read(const dictionary& dict);
private:
    fileName file_;

    typedef interpolationTable<scalar,polar*,2> polarTable;
    List<autoPtr<polar>> polars_;
    autoPtr<polarTable> polarInterpolated;

};


}


#endif /* PROPELLER_AIRFOIL_MODEL_POLAR_AIRFOIL_POLAR_AIRFOIL */
