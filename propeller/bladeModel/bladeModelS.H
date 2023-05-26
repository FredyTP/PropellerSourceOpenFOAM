#ifndef PROPELLER_BLADE_MODEL_BLADE_MODEL_S
#define PROPELLER_BLADE_MODEL_BLADE_MODEL_S

#include "linearInterpolation.H"
#include "airfoilModelList.H"
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

/**
 * @brief Class that represents a wing/blade, 
 * has 4 properties:
 *  - chord
 *  - twist
 *  - sweep
 *  - airfoil
 * The blade is made from diferent secctions or specifying the 
 * properties from a csv a list or constant value, linear or cubic spline
 * can be selected for the geometric properties (chord, twist, sweep).
 * 
 * The twist and sweep can be set from radian or degree, working internally in radian
 * 
 */
class bladeModelS
{
protected:
        //-------------------------------------//
        //------------BLADE PROPERTIES-----------//

        //- Airfoil of each section
        autoPtr<regularInterpolation<scalar,scalar,1>> chord_;
        autoPtr<regularInterpolation<scalar,scalar,1>> twistAngle_;
        autoPtr<regularInterpolation<scalar,scalar,1>> sweepAngle_;
        linearInterpolation<scalar,const airfoilModel*,1> airfoils_;

        bool isRadian_ = false;

        //- Optional input file
        fileName fName_;

        void readField(const dictionary& dict,csvTable<scalar,word>& table,autoPtr<regularInterpolation<scalar,scalar,1>>& field, bool isAngle = false);

        void checkBlade();
        void checkRadiusList(const List<scalar>& radiuslist);

        void readFromSections(const airfoilModelList &airfoilList, const dictionary& dict);
        void readFromProperties(const airfoilModelList &airfoilList, const dictionary& dict);
public:

    bladeModelS(const airfoilModelList& airfoilList, const dictionary& dict);


    bool sectionAtRadius(scalar adimRadius, scalar& chord, scalar& twist, scalar& sweep,interpolatedAirfoil& airfoil) const;
    bool geometryAtRadius(scalar adimRadius, scalar& chord, scalar& twist, scalar& sweep) const;

    /**
     * @brief Return unchecked sweep
     * 
     * @param adimRadius 
     * @return scalar 
     */
    scalar getSweep(scalar adimRadius) const;
    void writeBlade(label np,fileName path);

};


}




#endif /* PROPELLER_BLADE_MODEL_BLADE_MODEL_S */

