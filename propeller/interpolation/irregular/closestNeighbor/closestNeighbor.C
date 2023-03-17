#include "closestNeighbor.H"
#include "label.H"

namespace Foam
{

    template <class typeIn, class typeOu, label dim>
    closestNeighbor<typeIn, typeOu, dim>::closestNeighbor(const List<FixedList<typeIn, dim>> inputs_, const List<typeOu> outputs_)
    : inputs(inputs_), outputs(outputs_)
    {
        //initialize some internal process data if required (?)
        if(inputs.size()!=outputs.size())
        {
            //Error size must be the same
        }
    }

    template <class typeIn, class typeOu, label dim>
    interpolated<typeIn, typeOu> closestNeighbor<typeIn, typeOu, dim>::interpolate(

        FixedList<typeIn, dim> input) const
    {
        interpolated<typeIn, typeOu> result;

        scalar minDistanceSq = 1e300;
        label minIndex = -1;

        // TODO: extract this to a function of findClosest
        for (label i = 0; i < inputs.size(); ++i)
        {

            const FixedList<typeIn, dim> &testPoint =
                inputs[i];

            scalar sqDist = irregularInterpolation<typeIn, typeOu, dim>::SqrDistance(input, testPoint);

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
        result.points[0] = outputs[minIndex];

        result.coefficients().resize(1);
        result.coefficients()[0] = 1;
}


}

