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
            const word interpolation,
            List<scalar>& alpha,
            List<scalar>& cl,
            List<scalar>& cd,
            scalar Re,
            scalar Ma
        ),
        (interpolation, alpha,cl,cd, Re, Ma)
    );

    polar(const word interpolation, List<scalar>& alpha, List<scalar>& cl, List<scalar>& cd, scalar Re,scalar Ma);

    polar(const word interpolation, const fileName filename, scalar Re, scalar Ma);

    virtual scalar cl(scalar alpha) const;
    virtual scalar cd(scalar alpha) const;

    scalar reynolds() const;
    scalar mach() const;

    virtual bool valid();

    ~polar();
    
    static autoPtr<Foam::polar> New(const word modelType,const word interpolation, const fileName filename, scalar Re, scalar Ma);
    static autoPtr<Foam::polar> New(const word modelType,const word interpolation, List<scalar>& alpha, List<scalar>& cl, List<scalar>& cd, scalar Re,scalar Ma);
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
