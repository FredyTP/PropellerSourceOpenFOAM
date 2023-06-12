#ifndef PROPELLER_UTIL_READ_HELPER
#define PROPELLER_UTIL_READ_HELPER
#include "dictionary.H"
#include "scalar.H"
#include "word.H"
#include "csvTable.H"
#include "interpolationTable.H"

#include "autoPtr.H"

namespace Foam
{

namespace util
{

template<class typeIn, class typeOu>
void ReadInterpolationList
(
    const dictionary& dict,
    word x_name,
    word y_name,
    List<typeIn>& xList,
    List<typeOu>& yList,
    bool convertToRad = false,
    bool enableCSV=false,
    const csvTable<scalar,word>* csv = nullptr
);

template<class typeIn, class typeOu>
autoPtr<interpolationTable<typeIn,typeOu,1>> NewInterpolationFromDict
(
    const dictionary& dict,
    const List<typeIn>& xList,
    const List<typeOu>& yList
);

autoPtr<interpolationTable<scalar,scalar,1>> NewInterpolationFromDict
(
    const dictionary& dict,
    word x_name,
    word y_name,
    bool convertToRad = false,
    bool enableCSV=false,
    const csvTable<scalar,word>* csv = nullptr
);

}
}

#include "readHelper.C"

#endif /* PROPELLER_UTIL_READ_HELPER */
