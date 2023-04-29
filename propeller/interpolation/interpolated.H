#ifndef PROPELLER_INTERPOLATION_INTERPOLATED
#define PROPELLER_INTERPOLATION_INTERPOLATED
#include "List.H"
#include "scalar.H"
#include <functional>

namespace Foam
{


template<class typeIn,class typeOu>
class interpolated
{
protected:
    List<typeIn> coeff;
    List<typeOu> nodes;
public:
    interpolated(){}

    List<typeIn>& coefficients(){return coeff;}
    List<typeOu>& points(){return nodes;}

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
