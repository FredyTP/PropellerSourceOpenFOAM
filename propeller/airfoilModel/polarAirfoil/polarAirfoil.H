#ifndef PROPELLER_AIRFOIL_MODEL_POLAR_AIRFOIL_POLAR_AIRFOIL
#define PROPELLER_AIRFOIL_MODEL_POLAR_AIRFOIL_POLAR_AIRFOIL

#include "airfoilModel.H"
#include "scalar.H"
#include "dictionary.H"
#include "polar.H"
#include "Tuple2.H"
#include "interpolationTable.H"
namespace Foam
{


/**
 * @brief Airfoil made out of an 2D interpolation table of Reynolds and Mach
 * 
 */
: public airfoilModel
{
public:
    
    TypeName("polarAirfoil");

    polarAirfoil(word name, const dictionary& dict);

    virtual scalar cl(scalar alfaRad, scalar reynolds, scalar mach) const;
    virtual scalar cd(scalar alfaRad, scalar reynolds, scalar mach) const;

    bool read(const dictionary& dict);
    bool readFromPolars(bool cubicSpline, word extrapolation);

    /**
     * Requierements: 
     * - have columns named AoA, CL, CD, Re, Ma or specify names
     * - have all AoA sorted in ascending order
     * 
     * note:
     * - No need to have same AoA for every Mach / Re
     * - Data points dont have to be structured but it's recommended
    */
    bool readFromCSV(bool cubicSpline, word extrapolation);

protected:
    inline word aoaName() const  {return colNames[0];}
    inline word clName() const {return colNames[1];}
    inline word cdName() const {return colNames[2];}
    inline word reName() const {return colNames[3];}
    inline word maName() const {return colNames[4];}

    scalar convertInternalRe(scalar externalRe) const;
    scalar convertExternalRe(scalar internalRe) const;
private:
    fileName file_;
    //Use logarithmic reynolds internally to interpolation
    bool logRe_;
    typedef InterpolationTable<scalar,polar*,2> polarTable;
    List<autoPtr<polar>> polars_;
    autoPtr<polarTable> polarInterpolated;
    
    FixedList<word,5> colNames{"AoA","CL","CD","Re","Ma"};
    
    bool isRadian_ = false;

};


}


#endif /* PROPELLER_AIRFOIL_MODEL_POLAR_AIRFOIL_POLAR_AIRFOIL */
