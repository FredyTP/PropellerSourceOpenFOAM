#ifndef PROPELLER_PROPELLER_MODEL_BLADE_ELEMENT_MODEL_EFFECTS3D_POLAR_CORRECTION
#define PROPELLER_PROPELLER_MODEL_BLADE_ELEMENT_MODEL_EFFECTS3D_POLAR_CORRECTION
#include "mathematicalConstants.H"
#include "dictionary.H"
#include "Enum.H"

namespace Foam
{
namespace effects3D
{
class polarCorrection
{
public:
    enum class model
    {
        None,
        Snel,
        SnelPumping
    };

    polarCorrection(const dictionary& dict);
    void correct(scalar& cl3d, scalar cl2d, scalar alpha, scalar chord, scalar radius, scalar maxRadius,scalar angularVelocity, scalar Veff);
    static const Enum<model> modelNames_;

/**
 * @brief Snell rotational correction on cl only
 * Implemented from (MODELLING ROTATIONAL)
 * Snel, H. (ECN), Houwink, R. (NLR), and
 *  Bosschers, J. (NLR) ;
 * ‘Sectional prediction of lift coefficients on rotat-
 * ing wind turbine blades in stall’.
 * ECN-C--93-052, Petten, December 1994
 * 
 * @param cl3d cl with 3d effects as return value
 * @param cl2d cl2d from polars
 * @param alpha angle of attack 
 * @param chord chord of the blade at radius
 * @param radius radius of the points
 * @param maxRadius max radius of the blades
 */
static void Snel(scalar& cl3d, scalar cl2d, scalar alpha, scalar chord, scalar radius, scalar maxRadius);
/**
 * @brief Snell extended rotational correction on cl only, with local speed-ratio 
 * from centrifugal pumping model
 * 
 * Implemented from (MODELLING ROTATIONAL)
 * Snel, H. (ECN), Houwink, R. (NLR), and
 * Bosschers, J. (NLR) ;
 * ‘Sectional prediction of lift coefficients on rotat-
 * ing wind turbine blades in stall’.
 * ECN-C--93-052, Petten, December 1994
 * 
 * @param cl3d cl with 3d effects as return value
 * @param cl2d cl2d from polars
 * @param alpha angle of attack 
 * @param chord chord of the blade at radius
 * @param radius radius of the points
 * @param maxRadius max radius of the blades
 * @param angularVelocity blade angularVelocity
 * @param Veff airflow velocity
 */
static void SnelPumping(scalar& cl3d, scalar cl2d, scalar alpha, scalar chord, scalar radius, scalar maxRadius,scalar angularVelocity, scalar Veff);

private:
    model usedModel_;

    constexpr static scalar deg30inRad = 30.0/180.0*constant::mathematical::pi;
    constexpr static scalar deg50inRad = 50.0/180.0*constant::mathematical::pi;
    constexpr static scalar deg20inRad = 20.0/180.0*constant::mathematical::pi;

}; 
}
}

#endif /* PROPELLER_PROPELLER_MODEL_BLADE_ELEMENT_MODEL_EFFECTS3D_POLAR_CORRECTION */
