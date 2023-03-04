#ifndef PROPELLER_INTERPOLATION_INTERPOLATION_TABLE
#define PROPELLER_INTERPOLATION_INTERPOLATION_TABLE

#include "scalar.H"
#include "FixedList.H"

namespace Foam
{

template<unsigned int dimension>
class interpolationTable
{
public:
    interpolationTable(){}
    
    virtual scalar interpolate(FixedList<scalar,dimension> inputValue) = 0;
    //Definition for single input tables
    scalar interpolate(scalar input)
    {
        return interpolate(FixedList<scalar,dimension>(input));
    }
    
};

}

#endif /* PROPELLER_INTERPOLATION_INTERPOLATION_TABLE */
