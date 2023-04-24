#include "bladeModelS.H"
#include "unitConversion.H"
#include "interpolateXY.H"
#include "IFstream.H"

void Foam::bladeModelS::readField(const dictionary &dict, csvTable<scalar,word>& table, linearInterpolation<scalar,scalar,1>& field, bool isAngle)
{
    word from = dict.get<word>("from");

    scalarList radiusList;
    scalarList valuesList;

    if(from == "csv")
    {
        word dataCol = dict.get<word>("data");
        word radiusCol = dict.get<word>("radius");

        valuesList = table.col(dataCol);
        radiusList = table.col(radiusCol);
    }
    else if(from == "list")
    {
        valuesList = dict.get<scalarList>("data");
        radiusList = dict.get<scalarList>("radius");

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

    field.setData({radiusList},valuesList);
    field.setExtrapolationMode(extrapolationMode::emZero);

}

void Foam::bladeModelS::checkBlade()
{
    scalar minRadius=-1;
    scalar maxRadius=-1;

    //First radius
    minRadius = chord_.getInputs()[0][0];
    maxRadius = chord_.getInputs()[0].last();

    checkRadiusList(maxRadius,minRadius,chord_.getInputs()[0]);
    checkRadiusList(maxRadius,minRadius,twistAngle_.getInputs()[0]);
    checkRadiusList(maxRadius,minRadius,sweepAngle_.getInputs()[0]);
    checkRadiusList(maxRadius,minRadius,airfoils_.getInputs()[0]);

    if(!adimensional_)
    {
        maxRadius_ = maxRadius;
    }
    else
    {
        maxRadius_ = NO_RADIUS;
    }

}

void Foam::bladeModelS::checkRadiusList(scalar maxRadius, scalar minRadius, const scalarList &radiuslist)
{
    if(maxRadius != radiuslist.last())
    {
        FatalErrorInFunction
            <<"Max radius is not consistent:"
            <<endl
            <<radiuslist
            <<exit(FatalError);
    }
    else if(minRadius != radiuslist.first())
    {
        FatalErrorInFunction
            <<"Min radius is not consistent:"
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

Foam::bladeModelS::bladeModelS(
    const airfoilModelList &airfoilList,
    const dictionary &dict)
{

    /*List<Tuple2<word,FixedList<scalar,4>>> sections;
    
    List<bladeSection> bladeSections;
    List<scalar> radius;

    
    fName_=dict.getOrDefault<fileName>("file","");
    if(!fName_.empty())
    {
        Foam::IFstream is(fName_);
        is  >> sections;
    }
    else
    {
        dict.readEntry("sections",sections);       
    }*/
    Info<<endl;
    adimensional_ = dict.getOrDefault<bool>("adimensionalRadius",false);
    isRadian_ = dict.getOrDefault<bool>("isRadian",false);
    innerRadius_ = dict.getOrDefault<scalar>("innerRadius",0.0);
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
    airfoils_.setExtrapolationMode(extrapolationMode::emZero);

    Info<<"Chord: " <<chord_<<endl;
    Info<<"Twist: " <<twistAngle_<<endl;
    Info<<"Sweep: " <<sweepAngle_<<endl;

    checkBlade();

}

bool Foam::bladeModelS::sectionAtRadius(scalar radius, scalar& chord, scalar& twist, scalar& sweep,interpolatedAirfoil& airfoil)
{

    if(adimensional_)
    {
        radius /= maxRadius_;
    }

    if(radius>1.0 || radius < innerRadius_)
    {
        chord=0.0;
        return false;
    }
    
    chord = chord_.interpolate({radius}).value();
    twist = twistAngle_.interpolate({radius}).value();
    sweep = sweepAngle_.interpolate({radius}).value();
    interpolated<scalar,const airfoilModel*> afl = airfoils_.interpolate({radius});

    airfoil = interpolatedAirfoil(afl);

    return true;
}

void Foam::bladeModelS::setMaxRadius(scalar radius)
{
    maxRadius_=radius;
}

Foam::scalar Foam::bladeModelS::maxRadius() const
{
    return maxRadius_;
}

