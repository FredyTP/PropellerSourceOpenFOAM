#include "ClosestNeighbor.H"
#include "label.H"

namespace Foam
{

template <class typeIn, class typeOu, label dim>
ClosestNeighbor<typeIn, typeOu, dim>::ClosestNeighbor(const List<FixedList<typeIn, dim>> inputs_, const List<typeOu> outputs_)
: inputs(inputs_), outputs(outputs_)
{
    //initialize some internal process data if required (?)
    if(inputs.size()!=outputs.size())
    {
        //Error size must be the same
    }
}

template <class typeIn, class typeOu, label dim>
Interpolated<typeIn, typeOu> ClosestNeighbor<typeIn, typeOu, dim>::interpolate(

    FixedList<typeIn, dim> input) const
{
    Interpolated<typeIn, typeOu> result;

    scalar minDistanceSq = 1e300;
    label minIndex = -1;

    // TODO: extract this to a function of findClosest
    for (label i = 0; i < inputs.size(); ++i)
    {

        const FixedList<typeIn, dim> &testPoint =
            inputs[i];

        scalar sqDist = IrregularInterpolation<typeIn, typeOu, dim>::SqrDistance(input, testPoint);

        if (sqDist < minDistanceSq)
        {
            minDistanceSq = sqDist;
            minIndex = i;
        }
        else if (sqDist == 0.0)
        {
            minIndex = i;
            break;
        }
    }

    result.points().resize(1);
    result.points()[0] = outputs[minIndex];

    result.coefficients().resize(1);
    result.coefficients()[0] = 1;

    return result;
}

template <class typeIn, class typeOu, label dim>
label ClosestNeighbor<typeIn, typeOu, dim>::size()
{
    return outputs.size();
}
template <class typeIn, class typeOu, label dim>
bool ClosestNeighbor<typeIn, typeOu, dim>::setRawData(List<List<typeIn>>& inputs_, List<typeOu>& outputs_)
{
    return false;
}

}