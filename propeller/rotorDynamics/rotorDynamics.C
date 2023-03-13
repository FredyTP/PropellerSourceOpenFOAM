#include "rotorDynamics.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

    defineTypeNameAndDebug(rotorDynamics,0);

    //Define run time table for selecting derived types
    defineRunTimeSelectionTable(rotorDynamics, dictionary);   

rotorDynamics::rotorDynamics()
{

}

Foam::autoPtr<Foam::rotorDynamics> rotorDynamics::New(const dictionary &dict)
{
     //Get model Type name (Ex: simpleAirfoil) 
    //From typeNkey from dictionary (airfoilModel)
    const word modelType(dict.get<word>("type")); 

    Info<< " Selecting " << typeName << " " << modelType << endl;

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

    return autoPtr<Foam::rotorDynamics>(ctorPtr(dict));
}

}

