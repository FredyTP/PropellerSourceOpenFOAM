#ifndef PROPELLER_AIRFOIL_MODEL_POLAR_AIRFOIL_POLAR
#define PROPELLER_AIRFOIL_MODEL_POLAR_AIRFOIL_POLAR

#include <string>
#include "regularInterpolation.H"
#include "autoPtr.H"

namespace Foam
{

class polar
{
public:
    
    polar(const word interpolation, List<scalar> alpha, List<scalar> cl, List<scalar> cd, scalar Re,scalar Ma);

    polar(const word interpolation,fileName filename, scalar Re, scalar Ma);

    scalar cl(scalar alpha);
    scalar cd(scalar alpha);

    ~polar();
    
private:
    autoPtr<regularInterpolation<1>> cl_alpha;
    autoPtr<regularInterpolation<1>> cd_alpha;
    scalar reynolds_;
    scalar mach_;
};




}

#endif /* PROPELLER_AIRFOIL_MODEL_POLAR_AIRFOIL_POLAR */
