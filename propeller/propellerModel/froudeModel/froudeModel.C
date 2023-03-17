#include "froudeModel.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{
    //set and define the type name "typeName"
    defineTypeNameAndDebug(froudeModel,0);

    //Add to run time table to dynamically select the class based on
    //the dictionary propellerModel atribute
    addToRunTimeSelectionTable(propellerModel,froudeModel,dictionary);
}

Foam::froudeModel::froudeModel
(
    const dictionary& dict
) : propellerModel(dict,typeName)
{

    Info<<"Creating froude Model"<<endl;
}

 Foam::scalar Foam::froudeModel::radius() const
 {
    return 0;
 }
 