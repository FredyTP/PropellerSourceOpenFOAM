#ifndef PROPELLER_UTIL_OUTPUT_FIELDS_OUTPUT_FIELDS
#define PROPELLER_UTIL_OUTPUT_FIELDS_OUTPUT_FIELDS

#include "regIOobject.H"

namespace Foam
{

class outputFields
    : public regIOobject
{
public:
    enum class fields
    {
        aoa,
        inflowAngle,
        weights,
        selectedCells,
        sampledCells
    };




private:

}



}

#endif /* PROPELLER_UTIL_OUTPUT_FIELDS_OUTPUT_FIELDS */
