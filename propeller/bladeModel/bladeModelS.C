#include "bladeModelS.H"
#include "unitConversion.H"
#include "interpolateXY.H"
#include "IFstream.H"
#include "csvTable.H"
#include "OFstream.H"
#include "cubicSplineInterpolation.H"
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

void Foam::bladeModelS::readField(const dictionary &dict, csvTable<scalar, word> &table, autoPtr<regularInterpolation<scalar,scalar,1>> &field, bool isAngle)
{
    word from = dict.get<word>("from");

    scalarList radiusList;
    scalarList valuesList;
    bool isCubic=false;
    
    if(from == "csv")
    {
        word dataCol = dict.get<word>("data");
        word radiusCol = dict.get<word>("radius");
        valuesList = table.col(dataCol);
        radiusList = table.col(radiusCol);
        isCubic = dict.getOrDefault<bool>("cubicSpline",false);
    }
    else if(from == "list")
    {
        valuesList = dict.get<scalarList>("data");
        radiusList = dict.get<scalarList>("radius");
        isCubic = dict.getOrDefault<bool>("cubicSpline",false);
    }
    else if(from == "constant")
    {
        scalar data = dict.get<scalar>("data");
        valuesList.append(data);
        valuesList.append(data);

        radiusList.append(0.0);
        radiusList.append(1.0);
    }
    else
    {
        FatalErrorInFunction
            <<"From: "<<from<< " doesn't exit."
            <<exit(FatalError);
    }

    if(isAngle)
    {
        if(!isRadian_)
        {
            forAll(valuesList,i)
            {
                valuesList[i]=degToRad(valuesList[i]);
            }
        }
    }
    if(isCubic && radiusList.size() < 3)
    {
        Warning
        <<"Cannot use cubic spline if list.size() < 3 using linear interpolation on: "
        <<radiusList<<endl;
    }
    if(isCubic && radiusList.size()>2)
    {
        field = autoPtr<regularInterpolation<scalar,scalar,1>>::NewFrom<cubicSplineInterpolation>(radiusList,valuesList);
    }
    else
    {
        field = autoPtr<regularInterpolation<scalar,scalar,1>>::NewFrom<linearInterpolation<scalar,scalar,1>>(FixedList<scalarList,1>({radiusList}),valuesList);
    }

}

void Foam::bladeModelS::checkBlade()
{
    checkRadiusList(chord_->getInputs()[0]);
    checkRadiusList(twistAngle_->getInputs()[0]);
    checkRadiusList(sweepAngle_->getInputs()[0]);
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
            chord_ = autoPtr<regularInterpolation<scalar,scalar,1>>::NewFrom<cubicSplineInterpolation>(radius,chord);
            twistAngle_ = autoPtr<regularInterpolation<scalar,scalar,1>>::NewFrom<cubicSplineInterpolation>(radius,twist);
            sweepAngle_ = autoPtr<regularInterpolation<scalar,scalar,1>>::NewFrom<cubicSplineInterpolation>(radius,sweep);
        }
        else
        {
            chord_ = autoPtr<regularInterpolation<scalar,scalar,1>>::NewFrom<linearInterpolation<scalar,scalar,1>>(FixedList<scalarList,1>({radius}),chord);
            twistAngle_ = autoPtr<regularInterpolation<scalar,scalar,1>>::NewFrom<linearInterpolation<scalar,scalar,1>>(FixedList<scalarList,1>({radius}),twist);
            sweepAngle_ = autoPtr<regularInterpolation<scalar,scalar,1>>::NewFrom<linearInterpolation<scalar,scalar,1>>(FixedList<scalarList,1>({radius}),sweep);
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
    if(!csvFile.empty())
    {
        bladeCSV.readFile(csvFile);
    }

    const dictionary& chordDist = dict.subDict("chord");
    const dictionary& twistDist = dict.subDict("twist");
    const dictionary& sweepDist = dict.subDict("sweep");
    const dictionary& airfoilDist = dict.subDict("airfoils");

    readField(chordDist,bladeCSV,chord_);
    readField(twistDist,bladeCSV,twistAngle_,true);
    readField(sweepDist,bladeCSV,sweepAngle_,true);

    wordList airfoilNames = airfoilDist.get<wordList>("data");
    scalarList radius = airfoilDist.get<scalarList>("radius"); 

    List<const airfoilModel*> airfoilModels(airfoilNames.size());

    forAll(airfoilNames,i)
    {
        airfoilModels[i]=airfoilList.getAirfoil(airfoilNames[i]);
    }
    airfoils_.setData({radius},airfoilModels);
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
    interpolated<scalar,const airfoilModel*> afl = airfoils_.interpolate({adimRadius});
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
