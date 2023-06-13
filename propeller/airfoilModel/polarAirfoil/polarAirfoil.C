#include "polarAirfoil.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "LinearInterpolation.H"
#include "csvTable.H"

namespace Foam
{
    //set and define the type name "typeName"
    defineTypeNameAndDebug(polarAirfoil,0);

    //Add to run time table to dynamically select the class based on
    //the dictionary airfoilModel atribute

    addToRunTimeSelectionTable(airfoilModel,polarAirfoil,dictionary);

polarAirfoil::polarAirfoil(word name, const dictionary& dict)
:   airfoilModel(name,dict)
{
    if(!this->read(dict))
    {
        FatalIOErrorInFunction(dict)
        << "Could not read or invalid polar data "
        << exit(FatalIOError);
    }
}


bool polarAirfoil::read(const dictionary& dict)
{
    //Extrapolation mode: none, viterna, symmetry ... (?)
    word extrapolation = dict.getOrDefault<word>("extrapolation","polar");

    isRadian_ = dict.getOrDefault<bool>("isRadian",false);
    logRe_ = dict.getOrDefault<bool>("logReynolds",true);
    bool cubicSpline = dict.getOrDefault("cubicSpline",false);
    bool ok;
    ok = dict.readIfPresent("polars",file_);
    if(ok)
    {   
        ok= this->readFromPolars(cubicSpline, extrapolation);
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
        ok = this->readFromCSV(cubicSpline, extrapolation);
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
        ReMa[i][0]=convertInternalRe(polars_[i]->reynolds());
        ReMa[i][1]=polars_[i]->mach();
    }

    polarInterpolated = autoPtr<InterpolationTable<scalar,polar*,2>>::NewFrom<LinearInterpolation<scalar,polar*,2>>();
    polarInterpolated->setRawData(ReMa,polarList);

    return true;
}
bool polarAirfoil::readFromPolars(bool cubicSpline, word extrapolation)
{
    // Filename  - ( Re - Ma)
    List<Tuple2<word,FixedList<scalar,2>>> polarFiles;
    csvTable<scalar,word> csvReader(true);
    //Read airfoil data
    
    Foam::IFstream is(file_);
    is  >> polarFiles;

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
    
        auto ptrPolar = polar::New(extrapolation,cubicSpline,polarpath,Re,Ma,isRadian_);
        polars_[i].reset(ptrPolar.release());
    }

    return true;
}
bool polarAirfoil::readFromCSV(bool cubicSpline, word extrapolation)
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
        << exit(FatalError);
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
        auto ptrPolar = polar::New(extrapolation,cubicSpline,file_,0,0,isRadian_);
        polars_.resize(1);
        polars_[0].reset(ptrPolar.release());
        return true;
    }
    
    scalar maPolar = -1;//ma[0];
    scalar rePolar = -1;//re[0];

    List<List<scalar>> aoas;
    List<List<scalar>> cls;
    List<List<scalar>> cds;
    List<Tuple2<scalar,scalar>> reMa;

    label idx = 0;
    label maxSize = 0;
    for(label i = 0; i<aoa.size();i++)
    {   
        if(ma[i] != maPolar || re[i] != rePolar || i==0)
        {
            rePolar = re[i];
            maPolar = ma[i];
            label mareIdx = reMa.find({rePolar,maPolar});
            if(mareIdx==-1)
            {
                maxSize++;
                idx=maxSize-1;
                reMa.append({rePolar,maPolar});
                aoas.resize(maxSize);
                cls.resize(maxSize);
                cds.resize(maxSize);
            }
            else
            {
                idx=mareIdx;
            }
        }
        aoas[idx].append(aoa[i]);
        cls[idx].append(cl[i]);
        cds[idx].append(cd[i]);
    }

    for(label i = 0;i<aoas.size();i++)
    {
        //Create polar
        auto ptrPolar = polar::New(extrapolation,cubicSpline,aoas[i],cls[i],cds[i],reMa[i].first(),reMa[i].second(),isRadian_);
        if(ptrPolar->valid())
        {
            //Increment size by 1
            polars_.resize(polars_.size()+1);
            polars_[polars_.size()-1].reset(ptrPolar.release());

        }
        else
        {
            Info<<"PolarAirfoil: "<<this->airfoilName()
            <<" in Polar (Re = "<<reMa[i].first()<<", Ma = "<<reMa[i].second()
            <<") from file: "<<file_<< " - is invalid, descarting ..."<<endl;
        }
    }

    if(polars_.size()>0)
    {
        return true;
    }
    return false;
}
scalar polarAirfoil::convertInternalRe(scalar externalRe) const
{
    return logRe_ ? log(externalRe+VSMALL) : externalRe;
}
scalar polarAirfoil::convertExternalRe(scalar internalRe) const
{
    return logRe_ ? exp(internalRe) : internalRe;
}
scalar polarAirfoil::cl(scalar alfaRad, scalar reynolds, scalar mach) const
{
    return polarInterpolated->interpolate({convertInternalRe(reynolds),mach}).value<scalar>([=](scalar val,polar* p){return val*p->cl(alfaRad);});
}
scalar polarAirfoil::cd(scalar alfaRad, scalar reynolds, scalar mach) const
{
    return polarInterpolated->interpolate({convertInternalRe(reynolds),mach}).value<scalar>([=](scalar val,polar* p){return val*p->cd(alfaRad);});
}


}

