#include "bemControl.H"
#include "addToRunTimeSelectionTable.H"
namespace Foam
{
 
defineTypeNameAndDebug(bemControl,0);
addToRunTimeSelectionTable(rotorControl,bemControl, dictionary);

bemControl::bemControl(const dictionary& dict)
{
    
}



}