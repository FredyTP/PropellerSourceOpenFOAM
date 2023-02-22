#ifndef PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEARINTERPOLATION1
#define PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEARINTERPOLATION1

#include "List.H"
#include "scalar.H"
#include "regularInterpolation.H"
namespace Foam
{

class linearInterpolation1 
: public regularInterpolation<1>
{
protected:
    const List<scalar> inputs_;
    const List<scalar> outputs_;
public:

    //Constructors

        //- Construct from input and output list
        linearInterpolation1
        (
            const List<scalar> inputs,
            const List<scalar> outputs
        );
        
    virtual scalar interpolate(FixedList<scalar,1> input);

    

    /**
     * @return 0 if out of bounds 1 if succesful
    */
    label findIndex(scalar input,label& i1, label& i2) const;

    void interpolationCoefficients(scalar input, label i1, label i2, scalar& a1, scalar &a2) const;
       
};

}


#endif //PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEARINTERPOLATION1
