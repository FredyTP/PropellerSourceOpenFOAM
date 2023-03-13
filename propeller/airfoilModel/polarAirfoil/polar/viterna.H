#ifndef PROPELLER_AIRFOIL_MODEL_EXTRAPOLATION_VITERNA_VITERNA
#define PROPELLER_AIRFOIL_MODEL_EXTRAPOLATION_VITERNA_VITERNA


#include "polar.H"
#include "scalar.H"

namespace Foam
{

/**
 * Performs Viterna extrapolation with flat plat theory for large angles of attack
 * Implemented according to : The Effect of Flat Plate Theory Assumption in Post-Stall Lift
and Drag Coefficients Extrapolation with Viterna Method.
 * Author: Faisal Mahmuddin
*/
class viterna : public polar
{
public:
    TypeName("viterna");

    struct viterna_data
    {
        scalar alpha_stall_pos; //rad
        scalar alpha_stall_neg; //rad 
        scalar cl_stall_pos;
        scalar cl_stall_neg;
        scalar cd_stall_pos;
        scalar cd_stall_neg;
        scalar cd_max;

    };
    
    viterna(const word interpolation, List<scalar> alpha, List<scalar> cl, List<scalar> cd, scalar Re,scalar Ma);
    viterna(const word interpolation, const fileName filename, scalar Re, scalar Ma);

    virtual scalar cl(scalar alpha) const;
    virtual scalar cd(scalar alpha) const;


private:

    // - Coefficient functions
    scalar calcA1(scalar cd_max) const;
    scalar calcA2(scalar cd_max, scalar cl_stall, scalar alpha_stall) const;
    scalar calcB1(scalar cd_max) const;
    scalar calcB2(scalar cd_max, scalar cd_stall,scalar alpha_stall) const;


    // -Viterna extrapolation
    /**
     * Performs viterna extrapolation function for cl
     * @param alpha angle between alpha_stall - pi
    */
    scalar viterna_cl(scalar A1, scalar A2, scalar alpha) const;
    /**
     * Performs viterna extrapolation function for cd
     * @param alpha angle between alpha_stall - pi
    */
    scalar viterna_cd(scalar B1, scalar B2, scalar alpha) const;

    // -Flat plate extrapolation
    scalar flat_plate_cl(scalar cl_max_flat_plate,scalar alpha) const;
    scalar flat_plate_cd(scalar cd_max_flat_plate,scalar alpha) const;

    viterna_data data_;

    // Coeficientes de la interpolacion
    scalar A1_p,A2_p,B1_p,B2_p; //Positivos
    scalar A1_n,A2_n,B1_n,B2_n; //Negativos

    void setParameters();
};

}
#endif /* PROPELLER_AIRFOIL_MODEL_EXTRAPOLATION_VITERNA_VITERNA */
