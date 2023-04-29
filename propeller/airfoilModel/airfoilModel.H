#ifndef PROPELLER_AIRFOIL_MODEL_AIRFOIL_MODEL
#define PROPELLER_AIRFOIL_MODEL_AIRFOIL_MODEL

#include "dictionary.H"
#include "scalar.H"
#include "autoPtr.H"
#include "interpolated.H"
namespace Foam
{

class airfoilModel
{
public:
    TypeName("airfoilModel");
    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            airfoilModel,
            dictionary,
            (
                const word name,
                const dictionary& dict
            ),
            (name, dict)
        );

    airfoilModel();
    airfoilModel(const word name);
    
    virtual ~airfoilModel() = default;
    
    const word& airfoilName() const {return name_;}
    
    virtual scalar cl(scalar alfaRad, scalar reynolds, scalar mach) const = 0;
    virtual scalar cd(scalar alfaRad, scalar reynolds, scalar mach) const = 0;

    static Foam::autoPtr<Foam::airfoilModel> New(const word name, const dictionary& dict);

    void writeAirfoil(fileName path, scalar alfaBegin, scalar alfaEnd, label nAlfa, scalar Reyn, scalar Mach);
    //virtual scalar cm(scalar alfaRad, scalar reynolds, scalar mach);
protected:

    const word name_;

};



}

#endif /* PROPELLER_AIRFOIL_MODEL_AIRFOIL_MODEL */
