#include "fixedPower.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    //set and define the type name "typeName"
    defineTypeNameAndDebug(fixedPower,0);

    //Add to run time table to dynamically select the class based on
    //the dictionary airfoilModel atribute

    addToRunTimeSelectionTable(rotorDynamics,fixedPower,dictionary);

    fixedPower::fixedPower()
    {
    }

    fixedPower::fixedPower(const dictionary &dict)
    {


}


}
