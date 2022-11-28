#ifndef PROPELLER_ROTOR_MESH_ROTOR_GEOMETRY
#define PROPELLER_ROTOR_MESH_ROTOR_GEOMETRY
#include "scalar.H"
#include "vector.H"

#define NO_RADIUS -1.0

namespace Foam
{
    struct rotorGeometry
    {
        scalar radius {NO_RADIUS};
        vector center {0,0,0};
        vector direction {0,0,0};

        bool isValid()
        {
            return ( radius != NO_RADIUS ) &&
                   ( magSqr(direction) > SMALL);
        }
    };
}



#endif /* PROPELLER_ROTOR_MESH_ROTOR_GEOMETRY */
