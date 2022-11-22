#ifndef PROPELLER_AIRFOIL_MODEL_AIRFOIL_MODEL_LIST
#define PROPELLER_AIRFOIL_MODEL_AIRFOIL_MODEL_LIST

#include "airfoilModel.H"
#include "PtrList.H"

namespace Foam
{


class airfoilModelList
: PtrList<airfoilModel>
{
private:
    
public:
    airfoilModelList(const dictionary& dict);

    const airfoilModel* getAirfoil(const word name) const;
};

}


#endif /* PROPELLER_AIRFOIL_MODEL_AIRFOIL_MODEL_LIST */
