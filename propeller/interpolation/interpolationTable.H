#ifndef PROPELLER_INTERPOLATION_INTERPOLATION_TABLE
#define PROPELLER_INTERPOLATION_INTERPOLATION_TABLE

#include "scalar.H"
#include "FixedList.H"
#include "interpolated.H"

namespace Foam
{


enum extrapolationMode
{
    emConstant = 0,
    emZero,
    emExtrapolate, //not in use
};
/**
 * Base class for interpolation tables with no specified interpolation algorithm
 * 
 * The template parameters are general and can be set to any class that implements defined operators
 * The required operators depend on the implementation but in general should be defined
 * typeIn : addition/substraction (+/-), multiplication operator (*), comparation (<) and equality (==)
 * typeOu : to be able to return interpolation value: addition (+) and multiplication(*)
 * dim : must be a positive integer
*/
template<class typeIn, class typeOu, label dim>
class interpolationTable
{
protected:
    extrapolationMode exMode = extrapolationMode::emConstant;
public:

    void setExtrapolationMode(extrapolationMode mode) {exMode = mode;}

    
    interpolationTable(){}
    
    /**
     * Return interpolated value from input
    */
    virtual interpolated<typeIn,typeOu> interpolate(FixedList<typeIn,dim> input) const = 0;

    /**
     * @brief returns a list of interpolated coefficients and nodes
     * 
     * The reconstructed interpolated value is : coeff_1 * value_1 + coeff_2 * value_2 ...
    */
    //virtual void interpolationList(FixedList<typeIn,dim> input,List<typeIn>& coeff,List<typeOu>& values)=0;

    //Definition for single input tables
    virtual interpolated<typeIn,typeOu> interpolate(typeIn input) const 
    {
        return interpolate(FixedList<typeIn,dim>(input));
    }

    virtual void setRawData(List<List<typeIn>>& inputs_, List<typeOu>& outputs_) = 0;

    
    
};

}

#endif /* PROPELLER_INTERPOLATION_INTERPOLATION_TABLE */
