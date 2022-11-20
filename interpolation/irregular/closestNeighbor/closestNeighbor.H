#ifndef PROPELLER_INTERPOLATION_IRREGULAR_CLOSEST_NEIGHBOR_CLOSEST_NEIGHBOR
#define PROPELLER_INTERPOLATION_IRREGULAR_CLOSEST_NEIGHBOR_CLOSEST_NEIGHBOR
#include "irregularInterpolation.H"

namespace Foam
{
template<unsigned int dimension>
class closestNeighbor
: public irregularInterpolation<dimension>
{
    public:

    closestNeighbor(
        const List<FixedList<scalar,dimension>>& input,
        const List<scalar>& output
    );
    
    virtual scalar interpolate(FixedList<scalar,dimension> inputValue);

};

//Add template definitions
#include "closestNeighbor.C"

}

#endif /* PROPELLER_INTERPOLATION_IRREGULAR_CLOSEST_NEIGHBOR_CLOSEST_NEIGHBOR */
