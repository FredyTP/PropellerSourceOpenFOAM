#ifndef PROPELLER_AIRFOIL_MODEL_EXTRAPOLATION_EXTRAPOLATION
#define PROPELLER_AIRFOIL_MODEL_EXTRAPOLATION_EXTRAPOLATION
namespace Foam
{
/**
 * Interface to airfoil extrapolation methods
*/
class extrapolation
{
public:

    virtual scalar cl(scalar alpha) = 0;
    virtual scalar cd(scalar alpha) = 0;
}



}

#endif /* PROPELLER_AIRFOIL_MODEL_EXTRAPOLATION_EXTRAPOLATION */
