#include "airfoilModel.H"
#include "runTimeSelectionTables.H"
#include "csvTable.H"
#include "OFstream.H"

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

void Foam::airfoilModel::writeAirfoil(fileName path,scalar alfaBegin, scalar alfaEnd, label nAlfa, scalar Reyn, scalar Mach)
{
    scalar da = (alfaEnd-alfaBegin)/(nAlfa-1);
    csvTable<scalar,word> csv(true);
    List<word> title({"alphaRad","cl","cd"});
    csv.setHeader(title);
    for(label i = 0;i < nAlfa;i++)
    {
        scalar alfa = alfaBegin + i*da;
        scalar cl = this->cl(alfa,Reyn,Mach);
        scalar cd = this->cd(alfa,Reyn,Mach);

        List<scalar> row({alfa,cl,cd});
        csv.addRow(row);
    }

    OFstream file(path);
    file<<csv;

}
