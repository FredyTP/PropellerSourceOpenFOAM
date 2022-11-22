#include "bladeModelDev.H"
#include "unitConversion.H"
#include "interpolateXY.H"
Foam::devel::bladeModel::bladeModel
(
    const airfoilModelList& airfoilList,
    const dictionary& dict
)
:
    airfoils_(),
    radius_(),
    chord_(),
    torsionAngle_(),
    sweepAngle_()
{

    List<Tuple2<word,FixedList<scalar,4>>> sections;
    dict.readEntry("sections",sections);

    if(sections.size())
    {
        airfoils_.setSize(sections.size());
        radius_.setSize(sections.size());
        chord_.setSize(sections.size());
        torsionAngle_.setSize(sections.size());
        sweepAngle_.setSize(sections.size());

        forAll(sections,i)
        {
            //Check for airfoil exist
            airfoils_[i] = airfoilList.getAirfoil(sections[i].first());

            //Check radius is incresing  radius[i]>radius[i-1]!!
            radius_[i] = sections[i].second()[0];
            chord_[i] = sections[i].second()[1];
            torsionAngle_[i] = degToRad(sections[i].second()[2]);
            sweepAngle_[i] = degToRad(sections[i].second()[3]);
        }

    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "No blade data specified"
            << exit(FatalIOError);
    }



}

scalar Foam::dev::bladeModel::chordAtRadius(scalar radius) const
{
    //TODO: optimize those functions
    if(radius<radius_[0])
    {
        return 0;
    }
    else if(radius>radius_[radius_.size()])
    {
        return 0;
    }
    return interpolationXY<scalar>(radius,radius_,chord_);
}
scalar Foam::dev::bladeModel::torsionAtRadius(scalar radius) const
{
    if(radius<radius_[0])
    {
        return 0;
    }
    else if(radius>radius_[radius_.size()])
    {
        return 0;
    }
    return interpolationXY<scalar>(radius,radius_,torsionAngle_);
}
scalar Foam::dev::bladeModel::sweepAtRadius(scalar radius) const
{
    if(radius<radius_[0])
    {
        return 0;
    }
    else if(radius>radius_[radius_.size()])
    {
        return 0;
    }
    return interpolationXY<scalar>(radius,radius_,sweepAngle_);
}
