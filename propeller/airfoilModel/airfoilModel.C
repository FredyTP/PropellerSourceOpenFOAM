#include "airfoilModel.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

    defineTypeNameAndDebug(airfoilModel,0);

    //Define run time table for selecting derived types
    defineRunTimeSelectionTable(airfoilModel, dictionary);   

}
Foam::airfoilModel::airfoilModel()
{
    
}
Foam::airfoilModel::airfoilModel(const word name)
    : name_(name)
{

}
Foam::autoPtr<Foam::airfoilModel> Foam::airfoilModel::New
(
    const word name,
    const dictionary& dict
)
{

    //Get model Type name (Ex: simpleAirfoil) 
    //From typeNkey from dictionary (airfoilModel)
    const word modelType(dict.get<word>("type")); 

    Info<< "    - Loading " << modelType <<": "<<name<< endl;

    //Find class contructor in tables
    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            typeName,
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<Foam::airfoilModel>(ctorPtr(name, dict));
}
