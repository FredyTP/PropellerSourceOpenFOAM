#ifndef PROPELLER_AIRFOIL_MODEL_POLAR_AIRFOIL_POLAR
#define PROPELLER_AIRFOIL_MODEL_POLAR_AIRFOIL_POLAR

#include "csvTableReader.H"
#include "interpolationTable.H"

namespace foam
{

class polar
{
public:
    autoPtr<interpolationTable> cl_alfa;
    autoPtr<interpolationTable> cd_alfa;
    scalar reynolds_;
    scalar mach_;

    polar(String interpolation, List<scalar> alpha, List<scalar> cl, List<scalar> cd, scalar Re,scalar Ma);

    polar(String interpolation,std::string filename, scalar Re, scalar Ma);

private:

}




}

#endif /* PROPELLER_AIRFOIL_MODEL_POLAR_AIRFOIL_POLAR */
