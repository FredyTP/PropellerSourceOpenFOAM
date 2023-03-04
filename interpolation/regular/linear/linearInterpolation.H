#ifndef PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEARINTERPOLATION1
#define PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEARINTERPOLATION1

#include "List.H"
#include "scalar.H"
#include "regularInterpolation.H"
namespace Foam
{

class linearInterpolation1D 
: public regularInterpolation1D
{
protected:
    const List<scalar> inputs_;
    const List<scalar> outputs_;
public:

    //Constructors

        //- Construct from input and output list
        linearInterpolation1D
        (
            const List<scalar> inputs,
            const List<scalar> outputs
        );
        
    virtual scalar interpolate(FixedList<scalar,1> input);

    
protected:
    /**
     * @return 0 if out of bounds 1 if succesful
    */
    label findIndex(scalar input,label& i1, label& i2) const;

    void interpolationCoefficients(scalar input, label i1, label i2, scalar& a1, scalar &a2) const;
public:
   
    static void InterpolationCoefficients(scalar input,const List<scalar>& inputs, label i1, label i2, scalar& a1, scalar &a2);
    static label FindIndex(scalar input, const List<scalar>& inputs, label &i1, label &i2);
    
};




}


#endif //PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEARINTERPOLATION1
