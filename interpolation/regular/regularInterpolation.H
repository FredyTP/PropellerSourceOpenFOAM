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
    const List<FixedList<scalar,dimension>> input_;
    const List<scalar> output_;
public:

    //Constructors

        //- Construct from input and output list
        regularInterpolation
        (
            const List<FixedList<scalar,dimension>> input,
            const List<scalar> output
        );
        
    virtual scalar interpolate(FixedList<scalar,dimension> input)=0;

};

}



#include "regularInterpolation.C"

#endif /* PROPELLER_INTERPOLATION_REGULAR_REGULAR_INTERPOLATION */
