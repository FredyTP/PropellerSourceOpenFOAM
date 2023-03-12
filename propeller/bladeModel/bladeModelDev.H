#ifndef PROPELLER_BLADE_MODEL_BLADE_MODEL_DEV
#define PROPELLER_BLADE_MODEL_BLADE_MODEL_DEV

#include "linearInterpolation.H"
#include "airfoilModelList.H"
#include "bladeSection.H"
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
        linearInterpolation<scalar,bladeSection,1> sections_;

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


    interpolatedBladeSection sectionAtRadius(scalar radius);
    /*scalar chordAtRadius(scalar radius) const;
    scalar twistAtRadius(scalar radius) const;
    scalar sweepAtRadius(scalar radius) const;*/

    scalar maxRadius() const;

    void setMaxRadius(scalar radius);

};


}


}

#endif /* PROPELLER_BLADE_MODEL_BLADE_MODEL_DEV */
