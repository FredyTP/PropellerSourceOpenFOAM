#ifndef PROPELLER_INTERPOLATION_REGULAR_REGULAR_INTERPOLATION
#define PROPELLER_INTERPOLATION_REGULAR_REGULAR_INTERPOLATION

#include "List.H"
#include "scalar.H"
#include "interpolationTable.H"
namespace Foam
{

template<label dimension>
class regularInterpolation 
: public interpolationTable<dimension>
{
protected:
public:
        
    virtual scalar interpolate(FixedList<scalar,dimension> input)=0;
        //Definition for single input tables
    virtual scalar interpolate(scalar input)
    {
        return interpolate(FixedList<scalar,1>(input));
    }

};

class regularInterpolation1D : public interpolationTable<1>
{

};

class regularInterpolation2D : public interpolationTable<2>
{

};

class regularInterpolation3D : public interpolationTable<3>
{

};


}



#include "regularInterpolation.C"

#endif /* PROPELLER_INTERPOLATION_REGULAR_REGULAR_INTERPOLATION */
