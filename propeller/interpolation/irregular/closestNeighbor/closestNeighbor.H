#ifndef PROPELLER_INTERPOLATION_IRREGULAR_CLOSEST_NEIGHBOR_CLOSEST_NEIGHBOR
#define PROPELLER_INTERPOLATION_IRREGULAR_CLOSEST_NEIGHBOR_CLOSEST_NEIGHBOR
#include "irregularInterpolation.H"

namespace Foam
{
template<class typeIn, class typeOu, label dim>
class closestNeighbor
: public irregularInterpolation<typeIn,typeOu,dim>
{
private:
    const List<FixedList<typeIn,dim>> inputs;
    const List<typeOu> outputs;
public:

    closestNeighbor(
        const List<FixedList<typeIn,dim>> inputs_,
        const List<typeOu> outputs_
    );
    
    virtual interpolated<typeIn,typeOu> interpolate(FixedList<typeIn,dim> input) const;

    virtual label size();

    bool setRawData(List<List<typeIn>>& inputs_, List<typeOu>& outputs_);
   
};
}
//Add template definitions
#include "closestNeighbor.C"



#endif /* PROPELLER_INTERPOLATION_IRREGULAR_CLOSEST_NEIGHBOR_CLOSEST_NEIGHBOR */
