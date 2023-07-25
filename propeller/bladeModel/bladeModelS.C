#include "bladeModelS.H"
#include "unitConversion.H"
#include "interpolateXY.H"
#include "IFstream.H"
#include "csvTable.H"
#include "OFstream.H"
#include "cubicSplineInterpolation.H"
#include "readHelper.H"

void Foam::bladeModelS::writeBlade(label np, fileName path)
{
    scalar dr = 1.0/(np-1);
    List<word> header = {"r","c","twist","sweep"};
    List<scalar> bladeSection(4);
    csvTable<scalar,word> csv(true);
    csv.setHeader(header);
    interpolatedAirfoil airofil;
    Info<<"creating csv.."<<endl;
    for(int i = 0; i< np;i++)
    {
        scalar r = i*dr;
        scalar c,t,s;
        bladeSection[0]=r;

        if(this->sectionAtRadius(r,c,t,s,airofil))
        {
            bladeSection[1]=c;
            bladeSection[2]=t;
            bladeSection[3]=s;
        }
        else
        {
            bladeSection[1]=0;
            bladeSection[2]=0;
            bladeSection[3]=0;
        }

        csv.addRow(bladeSection);
    }

    OFstream file(path);
    file<<csv;

}

void Foam::bladeModelS::checkBlade()
{
    using namespace constant::mathematical;
    checkRadiusList(chord_->getInputs()[0]);
    checkRadiusList(twistAngle_->getInputs()[0]);
    checkRadiusList(sweepAngle_->getInputs()[0]);

    //Check sweep
    const List<scalar>& sweeps = sweepAngle_->getOutputs();
    forAll(sweeps,i)
    {
        if(sweeps[i]>pi/2 || sweeps[i]<pi/2)
        {
            FatalErrorInFunction
                <<"Sweep angle  should be between -pi/2 and pi/2"
                <<endl
                <<sweeps[i]
                <<exit(FatalError);
        }
    }
    checkRadiusList(airfoils_.getInputs()[0]);
}

void Foam::bladeModelS::checkRadiusList(const List<scalar>& radiuslist)
{
    if(radiuslist.first()<0.0)
    {
        FatalErrorInFunction
            <<"Radius cannont be negative, should be between 0.0 and 1.0"
            <<endl
            <<radiuslist
            <<exit(FatalError);
    }
    else if(radiuslist.last()>1.0)
    {
        FatalErrorInFunction
            <<"Radius cannont be greater than 1.0, should be between 0.0 and 1.0"
            <<endl
            <<radiuslist
            <<exit(FatalError);
    }
    else if(!isMonotonic<scalar>(radiuslist))
    {
        FatalErrorInFunction
            <<"Radius does not grow monotonically"
            <<endl
            <<radiuslist
            <<exit(FatalError);
    }
}

void Foam::bladeModelS::readFromSections(const airfoilModelList &airfoilList, const dictionary &dict)
{
    List<Tuple2<word,FixedList<scalar,4>>> sections;
    
    List<scalar> radius;
    List<scalar> chord;
    List<scalar> twist;
    List<scalar> sweep;
    List<const airfoilModel*> airfoils;


    dict.readEntry("sections",sections);       
    bool isCubic = dict.getOrDefault<bool>("cubicSpline",false);

    if(sections.size())
    {
        radius.resize(sections.size());
        chord.resize(sections.size());
        twist.resize(sections.size());
        sweep.resize(sections.size());
        airfoils.resize(sections.size());
        
        forAll(sections,i)
        {
            //Check for airfoil exist
            airfoils[i] = airfoilList.getAirfoil(sections[i].first());

            radius[i] = sections[i].second()[0];
            chord[i] = sections[i].second()[1];
            twist[i] = isRadian_ ? sections[i].second()[2] : degToRad(sections[i].second()[2]);
            sweep[i] = isRadian_ ? sections[i].second()[3] : degToRad(sections[i].second()[3]);
        }

        if(isCubic && radius.size() < 3)
        {
            Warning
            <<"Cannot use cubic spline if list.size() < 3 using linear interpolation on: "
            <<radius<<endl;
        }
        if(isCubic && radius.size()>2)
        {
            chord_ = autoPtr<RegularInterpolation<scalar,scalar,1>>::NewFrom<cubicSplineInterpolation>(radius,chord);
            twistAngle_ = autoPtr<RegularInterpolation<scalar,scalar,1>>::NewFrom<cubicSplineInterpolation>(radius,twist);
            sweepAngle_ = autoPtr<RegularInterpolation<scalar,scalar,1>>::NewFrom<cubicSplineInterpolation>(radius,sweep);
        }
        else
        {
            chord_ = autoPtr<RegularInterpolation<scalar,scalar,1>>::NewFrom<LinearInterpolation<scalar,scalar,1>>(FixedList<scalarList,1>({radius}),chord);
            twistAngle_ = autoPtr<RegularInterpolation<scalar,scalar,1>>::NewFrom<LinearInterpolation<scalar,scalar,1>>(FixedList<scalarList,1>({radius}),twist);
            sweepAngle_ = autoPtr<RegularInterpolation<scalar,scalar,1>>::NewFrom<LinearInterpolation<scalar,scalar,1>>(FixedList<scalarList,1>({radius}),sweep);
        }
        airfoils_.setData({radius},airfoils);
        
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "No blade data specified"
            << exit(FatalIOError);
    }
}

void Foam::bladeModelS::readFromProperties(const airfoilModelList &airfoilList, const dictionary &dict)
{
    fileName csvFile = dict.getOrDefault<fileName>("csv","");

    csvTable<scalar,word> bladeCSV(true);
    csvTable<scalar,word>* ptrCSV = nullptr;
    if(!csvFile.empty())
    {
        bladeCSV.readFile(csvFile);
        ptrCSV = &bladeCSV;
    }

    const dictionary& chordDict = dict.subDict("chord");
    const dictionary& twistDict = dict.subDict("twist");
    const dictionary& sweepDict = dict.subDict("sweep");
    const dictionary& airfoilDict = dict.subDict("airfoil");


    chord_ = util::NewInterpolationFromDict
    (
        chordDict,
        "radius",
        "chord",
        false,
        true,
        ptrCSV
    );
    twistAngle_ = util::NewInterpolationFromDict
    (
        twistDict,
        "radius",
        "twist",
        !isRadian_,
        true,
        ptrCSV
    );
    sweepAngle_ = util::NewInterpolationFromDict
    (
        sweepDict,
        "radius",
        "sweep",
        !isRadian_,
        true,
        ptrCSV
    );

    scalarList radius = airfoilDict.get<scalarList>("radius"); 
    wordList airfoilNames = airfoilDict.get<wordList>("airfoil");

    List<const airfoilModel*> airfoilModels(airfoilNames.size());

    forAll(airfoilNames,i)
    {
        airfoilModels[i]=airfoilList.getAirfoil(airfoilNames[i]);
    }
    airfoils_.setData({radius},airfoilModels);
}

Foam::bladeModelS::bladeModelS(const List<scalar> &rAdim, const List<scalar> &chord)
{
    checkRadiusList(rAdim);
    chord_ = autoPtr<RegularInterpolation<scalar,scalar,1>>
    ::NewFrom<LinearInterpolation<scalar,scalar,1>>(FixedList<scalarList,1>({rAdim}),chord);
            
    twistAngle_= autoPtr<RegularInterpolation<scalar,scalar,1>>
    ::NewFrom<LinearInterpolation<scalar,scalar,1>>(FixedList<scalarList,1>({rAdim}),chord);
            
    sweepAngle_= autoPtr<RegularInterpolation<scalar,scalar,1>>
    ::NewFrom<LinearInterpolation<scalar,scalar,1>>(FixedList<scalarList,1>({rAdim}),chord);
            
}

Foam::bladeModelS::bladeModelS(
    const airfoilModelList &airfoilList,
    const dictionary &dict)
{

    isRadian_ = dict.getOrDefault<bool>("isRadian",false);
    word mode = dict.get<word>("from");
    if(mode == "sections")
    {
        this->readFromSections(airfoilList,dict);
    }
    else if(mode == "properties")
    {
        this->readFromProperties(airfoilList,dict);
    }
    else
    {
        FatalErrorInFunction<<"from mode: "
        <<mode<<" does not exit. Valids modes are: (sections properties)"
        <<exit(FatalError);
    }
    Info<<endl;
    Info<<"Building blade from "<<mode<<endl;

    //Constant extrapolation to all span
    chord_->setExtrapolationMode(extrapolationMode::emConstant);
    twistAngle_->setExtrapolationMode(extrapolationMode::emConstant);
    sweepAngle_->setExtrapolationMode(extrapolationMode::emConstant);
    airfoils_.setExtrapolationMode(extrapolationMode::emConstant);

    checkBlade();
}

bool Foam::bladeModelS::sectionAtRadius(scalar adimRadius, scalar& chord, scalar& twist, scalar& sweep,interpolatedAirfoil& airfoil) const
{

    bool inside = geometryAtRadius(adimRadius,chord,twist,sweep);

    if(!inside)
    {
        return false;
    }
    Interpolated<scalar,const airfoilModel*> afl = airfoils_.interpolate({adimRadius});
    airfoil = interpolatedAirfoil(afl);

    return true;
}

bool Foam::bladeModelS::geometryAtRadius(scalar adimRadius, scalar &chord, scalar &twist, scalar &sweep) const
{
    if(adimRadius <0.0 || adimRadius > 1.0)
    {
        return false;
    }

    chord = chord_->interpolate({adimRadius}).value();
    twist = twistAngle_->interpolate({adimRadius}).value();
    sweep = sweepAngle_->interpolate({adimRadius}).value();

    return true;
}

Foam::scalar Foam::bladeModelS::getSweep(scalar adimRadius) const
{
    return sweepAngle_->interpolate({adimRadius}).value();
}