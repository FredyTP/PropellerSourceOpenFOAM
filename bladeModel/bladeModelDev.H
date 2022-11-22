#ifndef PROPELLER_BLADE_MODEL_BLADE_MODEL
#define PROPELLER_BLADE_MODEL_BLADE_MODEL


#include "airfoilModelList.H"

#include "scalar.H"
#include "label.H"
#include "List.H"
#include "vector.H"
#include "Tuple2.H"
namespace Foam
{

namespace devel
{

class bladeModel
{
private:

    List<const airfoilModel*> airfoils_;

    List<scalar> radius_;       // [m]
    List<scalar> chord_;        // [m]
    List<scalar> torsionAngle_; // [rads]
    List<scalar> sweepAngle_;   // [rads]


public:

    bladeModel(const airfoilModelList& airfoilList, const dictionary& dict);

    scalar chordAtRadius(scalar radius) const;
    scalar torsionAtRadius(scalar radius) const;
    scalar sweepAtRadius(scalar radius) const;

};


}


}

#endif /* PROPELLER_BLADE_MODEL_BLADE_MODEL */
