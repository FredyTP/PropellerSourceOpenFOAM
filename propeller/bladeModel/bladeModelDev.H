#ifndef PROPELLER_BLADE_MODEL_BLADE_MODEL_DEV
#define PROPELLER_BLADE_MODEL_BLADE_MODEL_DEV


#include "airfoilModelList.H"

#include "scalar.H"
#include "label.H"
#include "List.H"
#include "vector.H"
#include "Tuple2.H"
#include "fileName.H"
#include "rotorGeometry.H"

namespace Foam
{

namespace devel
{

class bladeModel
{
protected:
    // Protected Data


        //-------------------------------------//
        //------------BLADE SECTIONS-----------//

        //- Airfoil of each section
        List<const airfoilModel*> airfoils_;

        //- Radius where the section is located in [m]
        List<scalar> radius_;

        //- Chort of each section in [m]
        List<scalar> chord_;        

        //- Twisted angle of each section in [rad]
        List<scalar> twistAngle_;

        //- Sweep angle of each section in [rad]
        List<scalar> sweepAngle_;  

        //- If radius is adimensional
        bool adimensional_ = false;

        //- Holds max radius if dimensional and
        //- specifies max radius for adimensional data
        //- adimensional Radius = radius / maxRadius
        scalar maxRadius_ = NO_RADIUS;

        //- Optional input file
        fileName fName_;


public:

    bladeModel(const airfoilModelList& airfoilList, const dictionary& dict);

    scalar chordAtRadius(scalar radius) const;
    scalar twistAtRadius(scalar radius) const;
    scalar sweepAtRadius(scalar radius) const;

    scalar maxRadius() const;

    void setMaxRadius(scalar radius);

};


}


}

#endif /* PROPELLER_BLADE_MODEL_BLADE_MODEL_DEV */
