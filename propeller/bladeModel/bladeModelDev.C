#include "bladeModelDev.H"
#include "unitConversion.H"
#include "interpolateXY.H"
#include "IFstream.H"

Foam::devel::bladeModel::bladeModel
(
    const airfoilModelList& airfoilList,
    const dictionary& dict
)
{

    List<Tuple2<word,FixedList<scalar,4>>> sections;
    
    List<bladeSection> bladeSections;
    List<scalar> radius;

    adimensional_ = dict.getOrDefault<bool>("adimensionalRadius",false);
    fName_=dict.getOrDefault<fileName>("file","");
    if(!fName_.empty())
    {
        Info<<"Reading blade data from: " << fName_<<endl;
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
        bladeSections.setSize(sections.size());
        radius.setSize(sections.size());

        forAll(sections,i)
        {
            //Check for airfoil exist
            bladeSections[i].airfoil = airfoilList.getAirfoil(sections[i].first());

            //Check radius is incresing  radius[i]>radius[i-1]!!
            radius[i] = bladeSections[i].radius = sections[i].second()[0];

            bladeSections[i].chord = sections[i].second()[1];
            bladeSections[i].twistAngle = degToRad(sections[i].second()[2]);
            bladeSections[i].sweepAngle = degToRad(sections[i].second()[3]);
        }

        if(!adimensional_)
        {
            maxRadius_ = bladeSections[bladeSections.size()-1].radius;
        }
        else
        {
            maxRadius_ = NO_RADIUS;
        }
        sections_.setData({radius},bladeSections);
        sections_.setExtrapolationMode(extrapolationMode::emZero);
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "No blade data specified"
            << exit(FatalIOError);
    }
}

Foam::interpolatedBladeSection Foam::devel::bladeModel::sectionAtRadius(scalar radius) 
{
    if(adimensional_)
    {
        radius/=(maxRadius_+VSMALL);
    }
    auto sec = sections_.interpolate({radius});
    return interpolatedBladeSection(sec);
}

void Foam::devel::bladeModel::setMaxRadius(scalar radius)
{
    maxRadius_=radius;
}
/*
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
}*/

Foam::scalar Foam::devel::bladeModel::maxRadius() const
{
    return maxRadius_;
}

