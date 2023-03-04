#include "linearInterpolation2D.H"

namespace Foam
{
linearInterpolation2D::linearInterpolation2D
(
    const List<scalar> inputs1,
    const List<scalar> inputs2,
    const List<List<scalar>> outputs
) : inputs1_(inputs1),inputs2_(inputs2),outputs_(outputs)
{
    
    label sizein1,sizein2, sizeout1,sizeout2;

    sizein1 = inputs1_.size();
    sizein2 = inputs2_.size();
    sizeout1 = outputs_.size();
    sizeout2 = outputs_[0].size();

    if(sizein1 * sizein2 != sizeout1 * sizeout2)
    {
        //Error in formatting
    }
}

scalar linearInterpolation2D::interpolate(FixedList<scalar, 2> input)
{
    label i1,i2,j1,j2;
    scalar x1,x2,y1,y2;

    label r = findIndexes(input,i1,i2,j1,j2);

    x1=inputs1_[i1];
    x2=inputs1_[i2];
    y1=inputs2_[j1];
    y2=inputs2_[j2];

    linearInterpolation1D x_y1({x1,x2},{outputs_[i1][j1],outputs_[i2][j1]});
    linearInterpolation1D x_y2({x1,x2},{outputs_[i1][j2],outputs_[i2][j2]});

    linearInterpolation1D x_y({y1,y2},{x_y1.interpolate({input[0]}),x_y2.interpolate({input[0]})});

    return x_y.interpolate({input[1]});
}

label linearInterpolation2D::findIndexes(FixedList<scalar, 2> input, label &i1, label &i2, label &j1, label &j2)
{
    label r1 = linearInterpolation1D::FindIndex(input[0],inputs1_,i1,i2);
    label r2 = linearInterpolation1D::FindIndex(input[1],inputs2_,j1,j2);

    return r1+r2;
}


}