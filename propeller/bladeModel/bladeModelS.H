#ifndef PROPELLER_BLADE_MODEL_BLADE_MODEL_S
#define PROPELLER_BLADE_MODEL_BLADE_MODEL_S

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
#include "csvTable.H"
#include "interpolatedAirfoil.H"
#include "regularInterpolation.H"

namespace Foam
{


class bladeModelS
{
protected:
        //-------------------------------------//
        //------------BLADE SECTIONS-----------//

        //- Airfoil of each section
        autoPtr<regularInterpolation<scalar,scalar,1>> chord_;
        autoPtr<regularInterpolation<scalar,scalar,1>> twistAngle_;
        autoPtr<regularInterpolation<scalar,scalar,1>> sweepAngle_;
        linearInterpolation<scalar,const airfoilModel*,1> airfoils_;


        //- If radius is adimensional
        bool adimensional_ = false;
        bool isRadian_ = false;
        //- Holds max radius if dimensional and
        //- specifies max radius for adimensional data
        //- adimensional Radius = radius / maxRadius
        scalar maxRadius_ = NO_RADIUS;
        scalar innerRadius_ = 0.0;

        //- Optional input file
        fileName fName_;

        

        void readField(const dictionary& dict,csvTable<scalar,word>& table,autoPtr<regularInterpolation<scalar,scalar,1>>& field, bool isAngle = false);

        void checkBlade();
        void checkRadiusList(scalar maxRadius, scalar minRadius, const scalarList& radiuslist);
public:

    bladeModelS(const airfoilModelList& airfoilList, const dictionary& dict);


    bool sectionAtRadius(scalar radius, scalar& chord, scalar& twist, scalar& sweep,interpolatedAirfoil& airfoil);
    /*scalar chordAtRadius(scalar radius) const;
    scalar twistAtRadius(scalar radius) const;
    scalar sweepAtRadius(scalar radius) const;*/

    scalar maxRadius() const;

    void setMaxRadius(scalar radius);

    void writeBlade(label np,fileName path);

};


}




#endif /* PROPELLER_BLADE_MODEL_BLADE_MODEL_S */

