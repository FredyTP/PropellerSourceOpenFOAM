#ifndef PROPELLER_UTIL_READ_HELPER
#define PROPELLER_UTIL_READ_HELPER
#include "dictionary.H"
#include "scalar.H"
#include "word.H"
#include "csvTable.H"
#include "RegularInterpolation.H"

#include "autoPtr.H"

namespace Foam
{

namespace util
{

autoPtr<RegularInterpolation<scalar,scalar,1>> NewInterpolationFromDict
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

#endif /* PROPELLER_UTIL_READ_HELPER */
