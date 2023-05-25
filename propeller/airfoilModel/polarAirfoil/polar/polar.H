#ifndef PROPELLER_AIRFOIL_MODEL_POLAR_AIRFOIL_POLAR
#define PROPELLER_AIRFOIL_MODEL_POLAR_AIRFOIL_POLAR

#include <string>
#include "regularInterpolation.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "scalar.H"


namespace Foam
{

class polar
{
public:
    TypeName("polar");

    declareRunTimeSelectionTable
    (
        autoPtr,
        polar,
        dictionary,
        (
            bool cubicSpline,
            List<scalar>& alpha,
            List<scalar>& cl,
            List<scalar>& cd,
            scalar Re,
            scalar Ma,
            bool isRadian
        ),
        (cubicSpline, alpha,cl,cd, Re, Ma, isRadian)
    );

    polar(bool cubicSpline, List<scalar>& alpha, List<scalar>& cl, List<scalar>& cd, scalar Re,scalar Ma, bool isRadian);

    polar(bool cubicSpline, const fileName filename, scalar Re, scalar Ma, bool isRadian);

    virtual scalar cl(scalar alpha) const;
    virtual scalar cd(scalar alpha) const;

    scalar reynolds() const;
    scalar mach() const;

    virtual bool valid();

    virtual ~polar() = default;
    
    static autoPtr<Foam::polar> New(const word modelType,bool cubicSpline, const fileName filename, scalar Re, scalar Ma, bool isRadian);
    static autoPtr<Foam::polar> New(const word modelType,bool cubicSpline, List<scalar>& alpha, List<scalar>& cl, List<scalar>& cd, scalar Re,scalar Ma, bool isRadian);
protected:

    void processData(List<scalar>& alpha, List<scalar>& cl, List<scalar>& cd);

    scalar alpha_min;
    scalar alpha_max;
    scalar cl_alpha_min;
    scalar cl_alpha_max;
    scalar cd_alpha_max;
    scalar cd_alpha_min;
private:
    autoPtr<regularInterpolation<scalar,scalar,1>> cl_alpha;
    autoPtr<regularInterpolation<scalar,scalar,1>> cd_alpha;
    scalar reynolds_;
    scalar mach_;
};




}

#endif /* PROPELLER_AIRFOIL_MODEL_POLAR_AIRFOIL_POLAR */
