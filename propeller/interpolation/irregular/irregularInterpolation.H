#ifndef PROPELLER_INTERPOLATION_IRREGULAR_IRREGULAR_INTERPOLATION
#define PROPELLER_INTERPOLATION_IRREGULAR_IRREGULAR_INTERPOLATION



#include "List.H"
#include "scalar.H"
#include "label.H"
#include "interpolationTable.H"
namespace Foam
{

template<class typeIn, class typeOu, label dim>
class irregularInterpolation 
: public interpolationTable<typeIn,typeOu,dim>
{
public:
    virtual interpolated<typeIn,typeOu> interpolate(FixedList<typeIn,dim> input) const=0;

    virtual ~irregularInterpolation() = default;
protected:

    static typeIn SqrDistance(const FixedList<typeIn,dim>& p1, const FixedList<typeIn,dim>& p2) 
    {
        typeIn sqDist= 0;
        
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
