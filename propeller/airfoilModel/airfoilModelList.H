#ifndef PROPELLER_AIRFOIL_MODEL_AIRFOIL_MODEL_LIST
#define PROPELLER_AIRFOIL_MODEL_AIRFOIL_MODEL_LIST

#include "airfoilModel.H"
#include "PtrList.H"

namespace Foam
{


/**
 * @brief List of airfoil, auto import from dictionary, owning them in a PtrList.
 * 
 */
class airfoilModelList
: PtrList<airfoilModel>
{
private:
    
public:
    airfoilModelList(const dictionary& dict);

    const airfoilModel* getAirfoil(const word name) const;
    const airfoilModel* getAirfoil(label index) const;
};

}


#endif /* PROPELLER_AIRFOIL_MODEL_AIRFOIL_MODEL_LIST */
