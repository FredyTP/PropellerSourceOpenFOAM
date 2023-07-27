#ifndef PROPELLER_INTERPOLATION_INTERPOLATION_TABLE
#define PROPELLER_INTERPOLATION_INTERPOLATION_TABLE

#include "scalar.H"
#include "FixedList.H"
#include "Interpolated.H"

namespace Foam
{


enum class extrapolationMode
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
class InterpolationTable
{
protected:
    extrapolationMode exMode = extrapolationMode::emConstant;
public:

    typedef InterpolationTable<typeIn,typeOu,dim> baseTable;

    void setExtrapolationMode(extrapolationMode mode) {exMode = mode;}

    
    InterpolationTable(){}

    virtual ~InterpolationTable() = default;
    
    /**
     * Return interpolated value from input
    */
    virtual Interpolated<typeIn,typeOu> interpolate(const FixedList<typeIn,dim>& input) const = 0;

    /**
     * @brief returns a list of interpolated coefficients and nodes
     * 
     * The reconstructed interpolated value is : coeff_1 * value_1 + coeff_2 * value_2 ...
    */
    //virtual void interpolationList(FixedList<typeIn,dim> input,List<typeIn>& coeff,List<typeOu>& values)=0;

    //Definition for single input tables
    virtual Interpolated<typeIn,typeOu> interpolate(typeIn input) const 
    {
        return interpolate(FixedList<typeIn,dim>(input));
    }

    /**
     * Set interpolation data from a set of raw points (x1 y1 z1 ...) => (o1)
     *                                                 (x2 y2 z2 ...) => (o1)
     *                                                        .            .
     *                                                        .            .
     *                                                        .            .
     * If data set is not correct, return false
     * if could be fit, return true
    */                                                         
    virtual bool setRawData(const List<List<typeIn>>& inputs_, const List<typeOu>& outputs_) = 0;

    virtual label size() = 0;

    /**
     * @brief Checks if the input list grown monotonically
     * 
     * @tparam type 
     * @param checkList 
     * @return true 
     * @return false inverseDistance
     */
    
    
};

template<class type>
    static bool isMonotonic(const List<type>& checkList)
    {
        for(label i=0;i<checkList.size()-1;i++)
        {
            if(checkList[i]>checkList[i+1])
            {
                return false;
            }
        }
        return true;
    }




}
#endif /* PROPELLER_INTERPOLATION_INTERPOLATION_TABLE */
