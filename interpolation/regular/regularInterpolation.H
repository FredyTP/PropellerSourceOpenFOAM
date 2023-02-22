#ifndef PROPELLER_INTERPOLATION_REGULAR_REGULAR_INTERPOLATION
#define PROPELLER_INTERPOLATION_REGULAR_REGULAR_INTERPOLATION

#include "List.H"
#include "scalar.H"
#include "interpolationTable.H"
namespace Foam
{

template<unsigned int dimension>
class regularInterpolation 
: public interpolationTable<dimension>
{
protected:
public:

    //Constructors

        regularInterpolation();
        
    virtual scalar interpolate(FixedList<scalar,dimension> input)=0;
        //Definition for single input tables
    virtual scalar interpolate(scalar input)
    {
        return interpolate(FixedList<scalar,1>(input));
    }

};

}



#include "regularInterpolation.C"

#endif /* PROPELLER_INTERPOLATION_REGULAR_REGULAR_INTERPOLATION */
