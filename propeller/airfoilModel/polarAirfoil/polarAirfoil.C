#include "polarAirfoil.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "linearInterpolation.H"
#include "csvTable.H"

namespace Foam
{
    //set and define the type name "typeName"
    defineTypeNameAndDebug(polarAirfoil,0);

    //Add to run time table to dynamically select the class based on
    //the dictionary airfoilModel atribute

    addToRunTimeSelectionTable(airfoilModel,polarAirfoil,dictionary);

polarAirfoil::polarAirfoil(const word name, const dictionary& dict)
:   airfoilModel(name)
{
    if(!this->read(dict))
    {
        FatalIOErrorInFunction(dict)
        << "Could not read or invalid polar data "
        << exit(FatalIOError);
    }
}

/*
bool polarAirfoil::readTable(const dictionary& dict)
{
    Info<<"Reading polar airfoil data for: " << this->airfoilName() << endl;
    bool ok=true;
    ok &= dict.readEntry("file",file_);
    word extrapolation = dict.getOrDefault<word>("extrapolation","polar");
    //Read airfoil data
    if(!file_.empty())
    {
        Info<<"Reading polar data from: "<<file_<<endl;
        csvTable<scalar,word> csv(true);
        csv.readFile(file_);

        auto Re = csv.col2("Re");
        auto cl = csv.col("Cl");
        auto cd = csv.col("Cd");
        auto aoa = csv.col("AoA")
        
    }
}*/

bool polarAirfoil::read(const dictionary& dict)
{

    Info<<"Reading polar airfoil data for: " << this->airfoilName() << endl;

    //Extrapolation mode: none, viterna, symmetry ... (?)
    word extrapolation = dict.getOrDefault<word>("extrapolation","polar");

    bool ok;
    ok = dict.readIfPresent("polars",file_);
    if(ok)
    {   
        ok= this->readFromPolars(extrapolation);
    }
    else
    {
        dict.readEntry("csv",file_);
        List<word> columnNames;
        dict.readIfPresent("col",columnNames);
        //Change first N names used
        forAll(columnNames,i)
        {
            if(i<colNames.size())
            {
                colNames[i]=columnNames[i];
            }
        }
        ok = this->readFromCSV(extrapolation);
    }
    
    if(!ok)
    {
        return false;
    }
    //Create list to create interpolation table
    List<polar*> polarList(polars_.size());
    List<List<scalar>> ReMa(polars_.size());
    forAll(polars_,i)
    {
        polarList[i]=polars_[i].get();
        ReMa[i].setSize(2);
        ReMa[i][0]=polars_[i]->reynolds();
        ReMa[i][1]=polars_[i]->mach();
    }

    polarInterpolated = autoPtr<interpolationTable<scalar,polar*,2>>::NewFrom<linearInterpolation<scalar,polar*,2>>();
    polarInterpolated->setRawData(ReMa,polarList);

    return true;
}
bool polarAirfoil::readFromPolars(word extrapolation)
{
    // Filename  - ( Re - Ma)
    List<Tuple2<word,FixedList<scalar,2>>> polarFiles;
    csvTable<scalar,word> csvReader(true);
    //Read airfoil data
    if(!file_.empty())
    {
        Info<<"Reading polar data from: "<<file_<<endl;
        Foam::IFstream is(file_);
        is  >> polarFiles;
    }
    polars_.resize(polarFiles.size());
    //Build airfoil polars
    forAll(polarFiles,i)
    {
        auto polarf = polarFiles[i];
        fileName polarfile = polarf.first();
        scalar Re = polarf.second()[0];
        scalar Ma = polarf.second()[1];

        fileName polarpath = file_;
        polarpath.replace_name(polarfile);
    
        auto ptrPolar = polar::New(extrapolation,"lineal",polarpath,Re,Ma);
        polars_[i].reset(ptrPolar.release());
    }

    

    return true;
}
bool polarAirfoil::readFromCSV(word extrapolation)
{
    
    csvTable<scalar,word> csvReader(true);
    csvReader.readFile(file_);

    if(csvReader.nCol() == 0)
    {
        return false;
    }
    List<scalar> aoa,cl,cd,ma,re;
    aoa = csvReader.col(aoaName());
    cl = csvReader.col(clName());
    cd = csvReader.col(cdName());
    re = csvReader.col(reName());
    ma = csvReader.col(maName());

    label len = aoa.size();

    //If no valid data (empty or inconsistent)
    if(len == 0 || len != cl.size() || len != cd.size())
    {   
        FatalErrorInFunction
        << "Cannot find some " 
        << colNames 
        << " columns."
        << exit(FatalIOError);
    }

    //If only Reynold provided, set mach to 0s
    if(ma.size() == 0 && re.size() == len)
    {
        ma.resize(len,0);
    }
    //If only mach provided, set reynolds to 0s
    else if(re.size() == 0 && ma.size() == len)
    {
        re.resize(len,0);
    }
    //If neither mach or reynolds, set both to 0s and read as single polar
    else
    {
        auto ptrPolar = polar::New(extrapolation,"lineal",file_,0,0);
        polars_.resize(1);
        polars_[0].reset(ptrPolar.release());
        return true;
    }

    List<scalar> aoaPolar;
    List<scalar> clPolar;
    List<scalar> cdPolar;
    scalar maPolar = ma[0];
    scalar rePolar = re[0];

    for(label i = 0; i<aoa.size();i++)
    {
        if(i == aoa.size()-1)
        {
            aoaPolar.append(aoa[i]);
            clPolar.append(cl[i]);
            cdPolar.append(cd[i]);
        }
        if(ma[i] != maPolar || re[i] != rePolar || i == aoa.size()-1)
        {

            //Create polar
            auto ptrPolar = polar::New(extrapolation,"lineal",aoaPolar,clPolar,cdPolar,rePolar,maPolar);
            if(ptrPolar->valid())
            {
                //Increment size by 1
                polars_.resize(polars_.size()+1);
                polars_[polars_.size()-1].reset(ptrPolar.release());

            }
            else
            {
                Info<<"PolarAirfoil: "<<this->name_
                <<" in Polar (Re = "<<rePolar<<", Ma = "<<maPolar
                <<") from file: "<<file_<< " - is invalid, descarting ..."<<endl;
            }

            maPolar = ma[i];
            rePolar = re[i];

            aoaPolar.clear();
            clPolar.clear();
            cdPolar.clear();
        }
        aoaPolar.append(aoa[i]);
        clPolar.append(cl[i]);
        cdPolar.append(cd[i]);
    }

    if(polars_.size()>0)
    {
        return true;
    }
    return false;
}
scalar polarAirfoil::cl(scalar alfaRad, scalar reynolds, scalar mach) const
{
    return polarInterpolated->interpolate({reynolds,mach}).value([=](scalar val,polar* p){return val*p->cl(alfaRad);});
}
scalar polarAirfoil::cd(scalar alfaRad, scalar reynolds, scalar mach) const
{
    return polarInterpolated->interpolate({reynolds,mach}).value([=](scalar val,polar* p){return val*p->cd(alfaRad);});
}


}

