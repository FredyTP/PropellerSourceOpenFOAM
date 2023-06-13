#ifndef PROPELLER_INTERPOLATION_REGULAR_CUBIC_SPLINE_INTERPOLATION
#define PROPELLER_INTERPOLATION_REGULAR_CUBIC_SPLINE_INTERPOLATION

#include "List.H"
#include "scalar.H"
#include "RegularInterpolation.H"
#include "SquareMatrix.H"
#include "simpleMatrix.H"

namespace Foam
{

/**
 * @brief Implemented from 
 * https://medium.com/eatpredlove/natural-cubic-splines-implementation-with-python-edf68feb57aa 
 * equations
 * 
 */
class cubicSplineInterpolation : public RegularInterpolation<scalar,scalar,1>
{
private:
    typedef RegularInterpolation<scalar,scalar,1> parent;
    scalarList bi,ci,di;
public:

    cubicSplineInterpolation()
    {

    }
    
    

    cubicSplineInterpolation(List<scalar> inputs_, List<scalar> outputs_) 
        : parent(FixedList<List<scalar>,1>({inputs_}),outputs_)
    {
        const List<scalar>& x = parent::inputs[0];
        const List<scalar>& y = outputs;
        label n = x.size()-1;
        SquareMatrix<scalar> matrix(n+1,0);
        scalarField col(n+1,0);

        scalarField d(n);
        scalarField D(n);
        di.resize(n);
        bi.resize(n);
        for(label i = 0; i<n;i++)
        {
            d[i] = x[i+1]-x[i];
            D[i] = y[i+1]-y[i];
        }

        matrix(0,0)=1;
        matrix(n,n)=1;

        for(label i = 1; i<n;i++)
        {
            matrix(i,i)=2*(d[i]+d[i-1]);
            matrix(i,i-1)=d[i-1];
            matrix(i,i+1)=d[i];
            col[i]=3*(D[i]/d[i]-D[i-1]/d[i-1]);
        }
        simpleMatrix<scalar> solver(matrix,col);
        ci = solver.solve();

        for(label i = 0; i <n;i++)
        {
            di[i]=(ci[i+1]-ci[i])/(3*d[i]);
            bi[i]=(D[i]/d[i]-d[i]/3*(2*ci[i]+ci[i+1]));
        }
    }
    cubicSplineInterpolation(FixedList<List<scalar>,1> inputs_, List<scalar> outputs_)
        : cubicSplineInterpolation(inputs_[0],outputs_)
    {

    }
    Interpolated<scalar,scalar> interpolate(FixedList<scalar,1> input) const
    {
        label i1,i2;
        Interpolated<scalar,scalar> interp;
        label result = RegularInterpolation<scalar,scalar,1>::FindIndex(input[0],parent::inputs[0],i1,i2);

        if(result==0)
        {                
            switch (this->exMode)
            {
            case extrapolationMode::emConstant :
                interp.points().append(parent::outputs[i1]);
                interp.coefficients().append(1);
                break;
            case extrapolationMode::emZero :
                interp.points().append(0);
                interp.coefficients().append(0);
                break;
            default:
                interp.points().append(0);
                interp.coefficients().append(0);
                break;
            }
            
        }
        else
        {
            scalar dx = input[0]-parent::inputs[0][i1];
            interp.points().append(parent::outputs[i1]+bi[i1]*dx+ci[i1]*dx*dx+di[i1]*dx*dx*dx);
            interp.coefficients().append(1);
        }


        return interp;

    }
    virtual label size() {
        label s = 0;
        forAll(this->size_,i)
        {
            s+=this->size_[i];
        }
        return s;
    }
};


}



#endif /* PROPELLER_INTERPOLATION_REGULAR_CUBIC_SPLINE_INTERPOLATION */
