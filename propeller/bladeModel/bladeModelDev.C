#include "bladeModelDev.H"
#include "unitConversion.H"
#include "interpolateXY.H"
#include "IFstream.H"

Foam::devel::bladeModel::bladeModel
(
    const airfoilModelList& airfoilList,
    const dictionary& dict
)
:
    airfoils_(),
    radius_(),
    chord_(),
    twistAngle_(),
    sweepAngle_()
{

    List<Tuple2<word,FixedList<scalar,4>>> sections;

    adimensional_ = dict.getOrDefault<bool>("adimensional",false);
    fName_=dict.getOrDefault<fileName>("file","");
    if(!fName_.empty())
    {
        Info<<"Reading blade data from: "<<fName_<<endl;
        Foam::IFstream is(fName_);
        is  >> sections;
    }
    else
    {
        Info <<"Reading blade data from dict: " << endl;
        dict.readEntry("sections",sections);       
    }

    if(sections.size())
    {
        Info <<"   Building blade ... " << endl;
        airfoils_.setSize(sections.size());
        radius_.setSize(sections.size());
        chord_.setSize(sections.size());
        twistAngle_.setSize(sections.size());
        sweepAngle_.setSize(sections.size());

        forAll(sections,i)
        {
            //Check for airfoil exist
            airfoils_[i] = airfoilList.getAirfoil(sections[i].first());

            //Check radius is incresing  radius[i]>radius[i-1]!!
            radius_[i] = sections[i].second()[0];
            chord_[i] = sections[i].second()[1];
            twistAngle_[i] = degToRad(sections[i].second()[2]);
            sweepAngle_[i] = degToRad(sections[i].second()[3]);
        }

        if(!adimensional_)
        {
            maxRadius_ = radius_[radius_.size()-1];
        }
        else
        {
            maxRadius_ = NO_RADIUS;
        }

    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "No blade data specified"
            << exit(FatalIOError);
    }
}

void Foam::devel::bladeModel::setMaxRadius(scalar radius)
{
    if(!adimensional_)
    {
        maxRadius_=radius;
    }
}
Foam::scalar Foam::devel::bladeModel::chordAtRadius(Foam::scalar radius) const
{
    //TODO: optimize those functions
    if(adimensional_ && maxRadius_!=NO_RADIUS)
    {
        radius /= maxRadius_;
    }
    if(radius<radius_[0])
    {
        return 0;
    }
    else if(radius>radius_[radius_.size()-1])
    {
        return 0;
    }
    return Foam::interpolateXY<Foam::scalar>(radius,radius_,chord_);
}
Foam::scalar Foam::devel::bladeModel::twistAtRadius(Foam::scalar radius) const
{
    if(adimensional_ && maxRadius_!=NO_RADIUS)
    {
        radius /= maxRadius_;
    }
    if(radius<radius_[0])
    {
        return 0;
    }
    else if(radius>radius_[radius_.size()-1])
    {
        return 0;
    }
    return Foam::interpolateXY<Foam::scalar>(radius,radius_,twistAngle_);
}
Foam::scalar Foam::devel::bladeModel::sweepAtRadius(Foam::scalar radius) const
{
    if(adimensional_ && maxRadius_!=NO_RADIUS)
    {
        radius /= maxRadius_;
    }
    if(radius<radius_[0])
    {
        return 0;
    }
    else if(radius>radius_[radius_.size()-1])
    {
        return 0;
    }
    return Foam::interpolateXY<Foam::scalar>(radius,radius_,sweepAngle_);
}

Foam::scalar Foam::devel::bladeModel::maxRadius() const
{
    return maxRadius_;
}