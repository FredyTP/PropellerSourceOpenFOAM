#ifndef PROPELLER_AIRFOIL_MODEL_SIMPLE_AIRFOIL_SIMPLE_AIRFOIL
#define PROPELLER_AIRFOIL_MODEL_SIMPLE_AIRFOIL_SIMPLE_AIRFOIL

#include "airfoilModel.H"
#include "scalar.H"
#include "dictionary.H"

namespace Foam
{

/**
 * @brief Implementation of flat plate equations
 * 
 */
class simpleAirfoil
: public airfoilModel
{
public:
    
    TypeName("simpleAirfoil");


    simpleAirfoil(word name, scalar cl0,scalar cl_max,scalar cd_max,scalar cd0);
    simpleAirfoil(word name, const dictionary& dict);

    virtual scalar cl(scalar alfaRad, scalar reynolds, scalar mach) const;
    virtual scalar cd(scalar alfaRad, scalar reynolds, scalar mach) const;

    bool read(const dictionary& dict);
private:

    scalar cl0_;
    scalar cd0_;
    scalar cd_max_;
    scalar cl_max_;

};


}


#endif /* PROPELLER_AIRFOIL_MODEL_SIMPLE_AIRFOIL_SIMPLE_AIRFOIL */
