#ifndef PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEAR_INTERPOLATION3D
#define PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEAR_INTERPOLATION3D
#include "regularInterpolation.H"
#include "linearInterpolation.H"
#include "label.H"
namespace Foam
{

class linearInterpolation3D
: public regularInterpolation3D
{
    protected:
    const List<scalar> inputs1_;
    const List<scalar> inputs2_;
    const List<scalar> inputs3_;
    const List<List<List<scalar>>> outputs_;
public:

    //Constructors

        //- Construct from input and output list
        linearInterpolation3D
        (
            const List<scalar> inputs1,
            const List<scalar> inputs2,
            const List<scalar> inputs3,
            const List<List<List<scalar>>> outputs
        );

        //- Construct from raw data (inp 1, inp 2, inp 3 out 3)
        linearInterpolation3D
        (
            const List<FixedList<scalar,4>> data
        );

        
    virtual scalar interpolate(FixedList<scalar,3> input);
};
}
#endif /* PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEAR_INTERPOLATION3D */




