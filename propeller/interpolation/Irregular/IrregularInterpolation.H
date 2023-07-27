#ifndef PROPELLER_INTERPOLATION_IRREGULAR_IRREGULAR_INTERPOLATION
#define PROPELLER_INTERPOLATION_IRREGULAR_IRREGULAR_INTERPOLATION



#include "List.H"
#include "scalar.H"
#include "label.H"
#include "InterpolationTable.H"
namespace Foam
{

template<class typeIn, class typeOu, label dim>
class IrregularInterpolation 
: public InterpolationTable<typeIn,typeOu,dim>
{
public:
    virtual Interpolated<typeIn,typeOu> interpolate(const FixedList<typeIn,dim>& input) const=0;

    virtual ~IrregularInterpolation() = default;
protected:

    static typeIn SqrDistance(const FixedList<typeIn,dim>& p1, const FixedList<typeIn,dim>& p2) 
    {
        typeIn sqDist= Zero;
        
        for(label j = 0 ; j < dim; ++j)
        {
            auto d = p2[j]-p1[j];
            sqDist+= (d*d);
        }
        return sqDist;
    }

};

}



#endif /* PROPELLER_INTERPOLATION_IRREGULAR_IRREGULAR_INTERPOLATION */
