#include "regularInterpolation.H"

namespace Foam
{


template<unsigned int dimension>
regularInterpolation<dimension>::regularInterpolation
(
    const List<FixedList<scalar,dimension>> input,
    const List<scalar> output
) : 
input_(input),
output_(output)
{
    if(input_.size() != output_.size())
    {
        //Error size must be equal
    }
}


}