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

    typedef interpolationTable<typeIn,typeOu,dim> baseTable;

    void setExtrapolationMode(extrapolationMode mode) {exMode = mode;}

    
    interpolationTable(){}

    virtual ~interpolationTable() = default;
    
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

    /**
     * Set interpolation data from a set of raw points (x1 y1 z1 ...) => (o1)
     *                                                 (x2 y2 z2 ...) => (o1)
     *                                                        .            .
     *                                                        .            .
     *                                                        .            .
     * If data set is not correct, return false
     * if could be fit, return true
    */                                                         
    virtual bool setRawData(List<List<typeIn>>& inputs_, List<typeOu>& outputs_) = 0;

    virtual label size() = 0;

    //static autoPtr<interpolationTable> New();
    static autoPtr<interpolationTable> NewFromRaw(List<List<typeIn>>& inputs_, List<typeOu>& outputs_);

    /**
     * @brief Checks if the input list grown monotonically
     * 
     * @tparam type 
     * @param checkList 
     * @return true 
     * @return false 
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

template<class typeIn, class typeOu, label dim>
class constantInterpolation : public interpolationTable<typeIn,typeOu,dim>
{
    public:
    constantInterpolation(const typeOu& value) : value_(value) {}

    bool setRawData(List<List<typeIn>>& inputs_, List<typeOu>& outputs_) { return true;}
    interpolated<typeIn,typeOu> interpolate(FixedList<typeIn,dim> input) const
    {
        return interpolated<typeIn,typeOu>({1},{value_});
    }
    label size() { return 1; }
    private:
    typeOu value_;
};

}

#endif /* PROPELLER_INTERPOLATION_INTERPOLATION_TABLE */
