#include "bemDirectControl.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
 
defineTypeNameAndDebug(bemDirectControl,0);
addToRunTimeSelectionTable(base_class,bemDirectControl, dictionary);


bemDirectControl::bemDirectControl(const dictionary& dict);


}