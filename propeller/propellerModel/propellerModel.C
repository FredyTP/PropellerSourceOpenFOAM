#include "propellerModel.H"





namespace Foam
{

    defineTypeNameAndDebug(propellerModel,0);

    //Define run time table for selecting derived types
    defineRunTimeSelectionTable(propellerModel, dictionary);   
}

Foam::propellerModel::propellerModel
(
    const dictionary& dict,
    const word& name
)
{


}