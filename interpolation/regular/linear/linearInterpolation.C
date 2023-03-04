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
    return linearInterpolation1D::FindIndex(input,inputs_,i1,i2);
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

label linearInterpolation1D::FindIndex(scalar input, const List<scalar> &inputs, label &i1, label &i2)
{
    if(input<inputs[0])
    {
        //indexes are both set to closest index
        i1=0;
        i2=0;
        //out of bounds
        return 0;
    }
    for(int i = 0; i< inputs.size()-1; i++)
    {
        //If (luckily) there is coincidence
        if(input == inputs[i])
        {
            i1=i;
            i2=i;
            return 1;
        }

        //Between i and i+1
        if(inputs[i] < input && input< inputs[i+1])
        {
            i1=i;
            i2=i+1;
            return 1;
        }
    }

    if(input==inputs[inputs.size()-1])
    {
        i1=inputs.size()-1;
        i2=i1;
        return 1;
    }

    //Outside upper bounds
    i1=inputs.size()-1;
    i2=i1;
    return 0;
}

}