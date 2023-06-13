#include "inverseDistance.H"
#include "label.H"

namespace Foam
{

template <class typeIn, class typeOu, label dim>
inverseDistance<typeIn, typeOu, dim>::inverseDistance(const List<FixedList<typeIn, dim>> inputs_, const List<typeOu> outputs_,label nIntPoints)
: inputs(inputs_), outputs(outputs_), nInterpolationPoints(nIntPoints)
{
    //initialize some internal process data if required (?)
    if(inputs.size()!=outputs.size())
    {
        //Error size must be the same
    }
}
template <class typeIn, class typeOu, label dim>
interpolated<typeIn, typeOu> inverseDistance<typeIn, typeOu, dim>::interpolate(

    FixedList<typeIn, dim> input) const
{

    // THIS ALGORITHM IS NOT OPTIMIZED!! IT CAN BE VERY SLOW IN BIG DATA BASES (O(N^2))
    interpolated<typeIn, typeOu> result;

    List<Tuple2<scalar,label>> dists(inputs.size());
    List<label> indexes;

    
    scalar minDistanceSq = 1e300;
    label minIndex = -1;

    // TODO: extract this to a function of findClosest
    for (label i = 0; i < inputs.size(); ++i)
    {
        const FixedList<typeIn, dim> &testPoint =
            inputs[i];
        dists[i].first() = sqrt(irregularInterpolation<typeIn, typeOu, dim>::SqrDistance(input, testPoint));    
        dists[i].second() = i;
    }

    std::sort(dists.begin(),dists.end(),
        [](Tuple2<scalar,label> a, Tuple2<scalar,label> b) 
        {
            return a.first()<b.first();
        }
    );

    //IF ANY POINT IS COINCIDENT
    if(dists[0].first() < VSMALL)
    {
        result.points().resize(1);
        result.points()[0] = outputs[dists[0].second()];

        result.coefficients().resize(1);
        result.coefficients()[0] = 1;

        return result;
    }

    //set to max size if n = -1
    label np =  nInterpolationPoints == -1?dists.size():nInterpolationPoints;
    result.points().resize(np);
    result.coefficients().resize(np);

    scalar sum=0;
    for(label i = 0 ; i < np;i++)
    {
        sum+= 1.0/dists[i].first();
        result.points()[i]=outputs[dists[i].second()]; // set output index
    }

    for(label i = 0 ; i < np;i++)
    {
        result.coefficients()[i] = 1/(sum*dists[i].first());
        result.points()[i]=outputs[dists[i].second()]; // set output index
    }

    return result;

        
}

template <class typeIn, class typeOu, label dim>
bool inverseDistance<typeIn, typeOu, dim>::setRawData(List<List<typeIn>> &inputs_, List<typeOu> &outputs_)
{
    /**
     * Can try to fix ill-formed data, or not...
    */
    if(inputs_.size()!=outputs_.size())
    {
        return false;
    }
    inputs.resize(inputs_.size());
    outputs.resize(inputs_.size());

    for(label i = 0; i< inputs.size();i++)
    {
        for(label j =0;j<dim;j++)
        {
            if(inputs_.size()<dim)
            {
                return false;
            }

            inputs[i][j]=inputs_[i][j];
        }
        outputs[i]=outputs_[i];
    }
    return true;
}

template <class typeIn, class typeOu, label dim>
label inverseDistance<typeIn, typeOu, dim>::size()
{
    return outputs.size();
}
}

