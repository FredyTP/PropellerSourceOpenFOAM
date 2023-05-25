#include "airfoilModel.H"
#include "runTimeSelectionTables.H"
#include "csvTable.H"
#include "OFstream.H"
#include "constants.H"
#include "interpolatedAirfoil.H"
namespace Foam
{

    defineTypeNameAndDebug(airfoilModel,0);

    //Define run time table for selecting derived types
    defineRunTimeSelectionTable(airfoilModel, dictionary);   

}
Foam::airfoilModel::airfoilModel()
{
    
}

Foam::airfoilModel::airfoilModel(word name)
    : name_(name), export_(false)
{
    
}
Foam::airfoilModel::airfoilModel(word name, const dictionary& dict)
    : name_(name)
{
    export_ = dict.getOrDefault<bool>("export",false);
    path_ = dict.getOrDefault<fileName>("path","");
    aoaBegin_ = dict.getOrDefault<scalar>("aoaBegin", -constant::mathematical::pi);
    aoaEnd_ = dict.getOrDefault<scalar>("aoaBegin", constant::mathematical::pi);
    
    nAoa_ = dict.getOrDefault<scalar>("aoaBegin", 360); //1deg res
    Re_ = dict.getOrDefault<List<scalar>>("Reynolds", {0}); 
    Ma_ = dict.getOrDefault<List<scalar>>("Mach", {0}); 
    join_ = dict.getOrDefault<bool>("join", true); 
}
Foam::autoPtr<Foam::airfoilModel> Foam::airfoilModel::New
(
    word name,
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

void Foam::airfoilModel::exportAirfoil()
{
    if(export_)
    {
        fileName base_name = path_.size()>0 ? path_ : static_cast<fileName>(name_);
        Info << indent << "- Exporting airfoil: "<<name_<<endl;
        if(!join_)
        {
            for(label iRe = 0 ; iRe < Re_.size();iRe++)
            {
                for(label iMa = 0 ; iMa < Ma_.size(); iMa++)
                {
                    fileName final_name = base_name + "Re_"+std::to_string(Re_[iRe]) + "Ma_"+std::to_string(Ma_[iMa])+".csv";

                    this->writeAirfoil(final_name,aoaBegin_,aoaEnd_,nAoa_,{Re_[iRe]},{Ma_[iMa]});
                }
            }
        }
        else
        {
            this->writeAirfoil(base_name+".csv",aoaBegin_,aoaEnd_,nAoa_,Re_,Ma_);
        }

    }
}

void Foam::airfoilModel::writeAirfoil(fileName path, scalar alfaBegin, scalar alfaEnd, label nAlfa, const List<scalar> &Reyn, const List<scalar> &Mach)
{
    scalar da = (alfaEnd-alfaBegin)/(nAlfa-1);
    
    label n = Reyn.size()*Mach.size()*nAlfa;
    List<scalar> clList(n),cdList(n),aoaList(n),reList(n),maList(n);

    label cont =0;
    for(label iRe = 0 ; iRe < Reyn.size();iRe++)
    {
        for(label iMa = 0 ; iMa < Mach.size();iMa++)
        {
            for(label i = 0; i < nAlfa;i++)
            {
                scalar alfa = alfaBegin + i*da;
                aoaList[cont] = alfa;
                clList[cont] = this->cl(alfa,Reyn[iRe],Mach[iMa]);
                cdList[cont] = this->cd(alfa,Reyn[iRe],Mach[iMa]);
                reList[cont] = Reyn[iRe];
                maList[cont] = Mach[iMa];
                ++cont;
            }
        }
    }
    csvTable<scalar,word> csv(true);

    csv.addCol(aoaList,"aoaRad");
    csv.addCol(clList,"cl");
    csv.addCol(cdList,"cd");
    csv.addCol(reList,"Re");
    csv.addCol(maList,"Ma");


    OFstream file(path);
    file<<csv;

}
