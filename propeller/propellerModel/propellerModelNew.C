#include "propellerModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::propellerModel> Foam::propellerModel::New
(
    const dictionary& dict
)
{
    //Get model Type name (Ex: froudeModel) 
    //From type key from dictionary (propellerModel)
    const word modelType(dict.get<word>("type")); 

    Info<< "    Selecting " << typeName << " " << modelType << endl;

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

    return autoPtr<Foam::propellerModel>(ctorPtr(dict));
}


// ************************************************************************* //
