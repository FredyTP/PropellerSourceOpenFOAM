#include "interpolationTable.H"

namespace Foam
{

template <class typeIn, class typeOu, label dim>
inline autoPtr<interpolationTable> interpolationTable<typeIn, typeOu, dim>::NewFromRaw(List<List<typeIn>> &inputs_, List<typeOu> &outputs_)
{
    
    //Priority order regular -> irregular
    autoPtr<interpolationTable> ptrTable = autoPtr<interpolationTable>::NewFrom<linearInterpolation>();

    if(ptrTable->setRawData(inputs_,outputs_))
    {
        return ptrTable;
    }

    ptrTable = autoPtr<interpolationTable>::NewFrom<inverseDistance>();

    if(ptrTable->setRawData(inputs_,outputs_))
    {
        return ptrTable;
    }

    ptrTable.clear();
    
    return ptrTable;
}

}