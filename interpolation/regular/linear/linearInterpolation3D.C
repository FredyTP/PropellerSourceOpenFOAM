#include "linearInterpolation3D.H"

namespace Foam
{
linearInterpolation3D::linearInterpolation3D
(
    const List<scalar> inputs1,
    const List<scalar> inputs2,
    const List<scalar> inputs3,
    const List<List<List<scalar>>> outputs
) : inputs1_(inputs1),inputs2_(inputs2),inputs3_(inputs3),outputs_(outputs)
{
    
    label sizein1,sizein2,sizein3, sizeout1,sizeout2,sizeout3;

    sizein1 = inputs1_.size();
    sizein2 = inputs2_.size();
    sizein3 = inputs3_.size();
    sizeout1 = outputs_.size();
    sizeout2 = outputs_[0].size();
    sizeout3 = outputs_[0][0].size();

    if(sizein1 * sizein2 *sizeout3!= sizeout1 * sizeout2*sizeout3)
    {
        //Error in formatting
    }
}

scalar linearInterpolation3D::interpolate(FixedList<scalar, 3> input)
{
    label i0,i1,j0,j1,k0,k1;
    scalar x0,x1,y0,y1,z0,z1;
    scalar xd,yd,zd;
    scalar x,y,z;

    x=input[0];
    y=input[1];
    z=input[2];

    label r0 = linearInterpolation<scalar,scalar,1>::FindIndex(x,inputs1_,i0,i1);
    label r1 = linearInterpolation<scalar,scalar,1>::FindIndex(y,inputs2_,j0,j1);
    label r2 = linearInterpolation<scalar,scalar,1>::FindIndex(z,inputs3_,k0,k1);
    
    x0=inputs1_[i0];
    x1=inputs1_[i1];

    y0=inputs2_[j0];
    y1=inputs2_[j1];

    z0=inputs3_[k0];
    z1=inputs3_[k1];

    xd=(x-x0)/(x1-x0);
    yd=(y-y0)/(y1-y0);
    zd=(x-x0)/(z1-z0);

    scalar o_y0z0 = outputs_[i0][j0][k0]*(1-xd) + outputs_[i1][j0][k0]*xd;
    scalar o_y1z0 = outputs_[i0][j1][k0]*(1-xd) + outputs_[i1][j1][k0]*xd;
    scalar o_y0z1 = outputs_[i0][j0][k1]*(1-xd) + outputs_[i1][j0][k1]*xd;
    scalar o_y1z1 = outputs_[i0][j1][k1]*(1-xd) + outputs_[i1][j1][k1]*xd;


    scalar o_z0 = o_y0z0*(1-yd) + o_y1z0*yd;
    scalar o_z1 = o_y0z1*(1-yd) + o_y1z1*yd;

    return o_z0*(1-zd)+o_z1*zd;
}




}