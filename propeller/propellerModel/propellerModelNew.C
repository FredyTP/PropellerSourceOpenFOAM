#include "propellerModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::propellerModel> Foam::propellerModel::New
(
    word sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)
{
    //Get model Type name (Ex: froudeModel) 
    //From type key from dictionary (propellerModel)
    const word modelType(dict.get<word>("type")); 

    Info<< "    -Selecting " << typeName << " " << modelType << endl;

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

    return autoPtr<Foam::propellerModel>(ctorPtr(sourceName,dict,mesh));
}


// ************************************************************************* //
