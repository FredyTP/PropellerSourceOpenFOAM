#ifndef PROPELLER_INTERPOLATION_IRREGULAR_CLOSEST_NEIGHBOR_CLOSEST_NEIGHBOR
#define PROPELLER_INTERPOLATION_IRREGULAR_CLOSEST_NEIGHBOR_CLOSEST_NEIGHBOR
#include "IrregularInterpolation.H"

namespace Foam
{
template<class typeIn, class typeOu, label dim>
class ClosestNeighbor
: public IrregularInterpolation<typeIn,typeOu,dim>
{
private:
    const List<FixedList<typeIn,dim>> inputs;
    const List<typeOu> outputs;
public:

    ClosestNeighbor(
        const List<FixedList<typeIn,dim>> inputs_,
        const List<typeOu> outputs_
    );
    
    virtual Interpolated<typeIn,typeOu> interpolate(FixedList<typeIn,dim> input) const;

    virtual label size();

    bool setRawData(List<List<typeIn>>& inputs_, List<typeOu>& outputs_);
   
};
}
//Add template definitions
#include "ClosestNeighbor.C"



#endif /* PROPELLER_INTERPOLATION_IRREGULAR_CLOSEST_NEIGHBOR_CLOSEST_NEIGHBOR */
