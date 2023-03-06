#include "linearInterpolation1D.H"
#include "linearInterpolation.H"

namespace Foam
{



linearInterpolation1D::linearInterpolation1D
(
    const List<scalar> inputs,
    const List<scalar> outputs
) : 
inputs_(inputs),
outputs_(outputs)
{
    if(inputs_.size() != outputs_.size())
    {
        //Error size must be equal
    }
}

scalar linearInterpolation1D::interpolate(FixedList<scalar, 1> input){
    scalar x = input[0];

    label i1,i2;

    label result = findIndex(x,i1,i2);

    //Out of bounds (extrapolate)
    if(result == 0)
    {
        //extrapolate
        return outputs_[i1];
    }
    //Coincident
    else if(i1==i2)
    {
        //
        return outputs_[i1];
    }
    else
    {
        scalar a1,a2;
        interpolationCoefficients(x,i1,i2,a1,a2);
        return a1*outputs_[i1] + a2 * outputs_[i2];
    }

}

label linearInterpolation1D::findIndex(scalar input,label& i1, label& i2) const
{
    return linearInterpolation<scalar,scalar,1>::FindIndex(input,inputs_,i1,i2);
}
void linearInterpolation1D::interpolationCoefficients(scalar input, label i1, label i2, scalar &a1, scalar &a2) const
{
    scalar x,x1,x2;

    x = input;
    x1 = inputs_[i1];
    x2 = inputs_[i2];

    a1 = (x2-x)/(x2-x1);
    a2 = (x-x1)/(x2-x1);
}

void linearInterpolation1D::InterpolationCoefficients(scalar input, const List<scalar> &inputs, label i1, label i2, scalar &a1, scalar &a2)
{
    scalar x,x1,x2;

    x = input;
    x1 = inputs[i1];
    x2 = inputs[i2];

    a1 = (x2-x)/(x2-x1);
    a2 = (x-x1)/(x2-x1);
}




}