#ifndef PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEAR_INTERPOLATION2D
#define PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEAR_INTERPOLATION2D
#include "regularInterpolation.H"
#include "linearInterpolation.H"
#include "label.H"
namespace Foam
{

class linearInterpolation2D
: public regularInterpolation2D
{
    protected:
    const List<scalar> inputs1_;
    const List<scalar> inputs2_;
    const List<List<scalar>> outputs_;
public:

    //Constructors

        //- Construct from input and output list
        linearInterpolation2D
        (
            const List<scalar> inputs1,
            const List<scalar> inputs2,
            const List<List<scalar>> outputs
        );

        //- Construct from raw data (inp 1, inp 2, out 3)
        linearInterpolation2D
        (
            const List<FixedList<scalar,3>> data
        );

        
    virtual scalar interpolate(FixedList<scalar,2> input);

    /**
     * @return 0 if out of bounds 1 if succesful
    */
    label findIndexes(FixedList<scalar,2> input,label& i1, label& i2, label& j1, label& j2);

    //void interpolationCoefficients(scalar input, label i1, label i2, scalar& a1, scalar &a2) const;

};
}
#endif /* PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEAR_INTERPOLATION2D */
