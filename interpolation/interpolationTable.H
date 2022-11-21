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
};


}

#endif /* PROPELLER_INTERPOLATION_INTERPOLATION_TABLE */
