#include "InterpolationTable.H"

namespace Foam
{

template <class typeIn, class typeOu, label dim>
inline autoPtr<InterpolationTable> InterpolationTable<typeIn, typeOu, dim>::NewFromRaw(List<List<typeIn>> &inputs_, List<typeOu> &outputs_)
{
    
    //Priority order regular -> irregular
    autoPtr<InterpolationTable> ptrTable = autoPtr<InterpolationTable>::NewFrom<linearInterpolation>();

    if(ptrTable->setRawData(inputs_,outputs_))
    {
        return ptrTable;
    }

    ptrTable = autoPtr<InterpolationTable>::NewFrom<inverseDistance>();

    if(ptrTable->setRawData(inputs_,outputs_))
    {
        return ptrTable;
    }

    ptrTable.clear();
    
    return ptrTable;
}

}