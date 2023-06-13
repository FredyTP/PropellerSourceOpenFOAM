#ifndef PROPELLER_INTERPOLATION_IRREGULAR_INVERSE_DISTANCE_INVERSE_DISTANCE
#define PROPELLER_INTERPOLATION_IRREGULAR_INVERSE_DISTANCE_INVERSE_DISTANCE

#include "IrregularInterpolation.H"

namespace Foam
{
/**
 * Inverse distance weighting algorithm for interpolating irregular data
*/
template<class typeIn, class typeOu, label dim>
class InverseDistance
: public IrregularInterpolation<typeIn,typeOu,dim>
{
private:
    const List<FixedList<typeIn,dim>> inputs;
    const List<typeOu> outputs;
    label nInterpolationPoints;
public:

    InverseDistance(
        const List<FixedList<typeIn,dim>> inputs_,
        const List<typeOu> outputs_,
        label nIntPoints = -1
    );
    
    virtual interpolated<typeIn,typeOu> interpolate(FixedList<typeIn,dim> input) const;

    virtual bool setRawData(List<List<typeIn>>& inputs_, List<typeOu>& outputs_);

    virtual label size();
   
};
}
//Add template definitions
#include "InverseDistance.C"



#endif /* PROPELLER_INTERPOLATION_IRREGULAR_INVERSE_DISTANCE_INVERSE_DISTANCE */

