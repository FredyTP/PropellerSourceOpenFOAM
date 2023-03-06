#include "closestNeighbor.H"



template<unsigned int dimension>
Foam::closestNeighbor<dimension>::closestNeighbor(
    const Foam::List<FixedList<scalar,dimension>> input,
    const Foam::List<scalar> output
)
: irregularInterpolation<dimension>(input,output)
{
    //initialize some internal process data if required (?)
}

template<unsigned int dimension>
Foam::scalar Foam::closestNeighbor<dimension>::interpolate
(

    Foam::FixedList<scalar,dimension> inputValue
)
{
    Foam::scalar minDistanceSq = 1e300;
    Foam::label minIndex=-1;

    //TODO: extract this to a function of findClosest
    forAll(Foam::irregularInterpolation<dimension>::input_,i)
    {
        Foam::scalar sqDist=0;
        const Foam::FixedList<scalar,dimension>& testPoint = 
        Foam::irregularInterpolation<dimension>::input_[i];

        //TODO: extract this to a function of compute N-distance
        forAll(inputValue,j)
        {
            sqDist+=Foam::sqr(inputValue[j]-testPoint[j]);
        }


        if(sqDist < minDistanceSq)
        {
            minDistanceSq=sqDist;
            minIndex=i;
        }
        else if(sqDist == 0.0)
        {
            minIndex = i;
            break;
        }
    }

    return Foam::irregularInterpolation<dimension>::output_[minIndex];
}

template<unsigned int dimension>
Foam::scalar Foam::closestNeighbor<dimension>::Interpolate
(
    const List<FixedList<scalar,dimension>>& inputList,
    const List<scalar>& outputList,
    FixedList<scalar,dimension> inputValue
)
{
    return 0;
}