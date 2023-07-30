#ifndef PROPELLER_INTERPOLATION_INTERPOLATED
#define PROPELLER_INTERPOLATION_INTERPOLATED
#include "List.H"
#include "scalar.H"
#include <functional>

namespace Foam
{

/**
 * @brief Class generalizing interpolation concept, enabled to contain list of values and coeffs
 * 
 * @tparam typeIn type of coefficientes
 * @tparam typeOu type of interpolated
 */
template<class typeIn,class typeOu>
class Interpolated
{
protected:
    List<typeIn> coeff;
    List<typeOu> nodes;
public:
    Interpolated(){}
    Interpolated(const List<typeIn>& coe, const List<typeOu>& nod) : coeff(coe), nodes(nod) {}
    List<typeIn>& coefficients(){return coeff;}
    List<typeOu>& points(){return nodes;}

    /**
     * @brief Templated method to obtaing interpolated valued
     * 
     * @tparam typeRet 
     * @param aggregate function that contains the logic to aggregate coeficients and values
     * @return typeRet 
     */
    template<class typeRet = typeOu>
    typeRet value(std::function<typeRet(typeIn,typeOu)> aggregate=
    [](typeIn a, typeOu b){return a*b;}) const
    {
        typeRet acc = Zero;
        for(label i = 0; i<coeff.size();i++)
        {
            acc += aggregate(coeff[i],nodes[i]);
        }
        return acc;
    }


};
}
#endif /* PROPELLER_INTERPOLATION_INTERPOLATED */
