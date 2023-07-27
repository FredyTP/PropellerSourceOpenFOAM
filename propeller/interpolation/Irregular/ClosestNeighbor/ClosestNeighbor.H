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
    List<FixedList<typeIn,dim>> inputs;
    List<typeOu> outputs;
public:

    ClosestNeighbor() = default;

    ClosestNeighbor(
        const List<FixedList<typeIn,dim>> inputs_,
        const List<typeOu> outputs_
    );
    
    virtual Interpolated<typeIn,typeOu> interpolate(const FixedList<typeIn, dim>& input) const;

    virtual label size();

    virtual bool setRawData(const List<List<typeIn>>& inputs_, const List<typeOu>& outputs_);
   
};
}
//Add template definitions
#include "ClosestNeighbor.C"



#endif /* PROPELLER_INTERPOLATION_IRREGULAR_CLOSEST_NEIGHBOR_CLOSEST_NEIGHBOR */
