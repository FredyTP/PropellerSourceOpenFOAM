#ifndef PROPELLER_AIRFOIL_MODEL_AIRFOIL_MODEL
#define PROPELLER_AIRFOIL_MODEL_AIRFOIL_MODEL

#include "dictionary.H"
#include "scalar.H"
#include "autoPtr.H"

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

    airfoilModel(const word name);
    
    
    const word& airfoilName(){return name_;}
    
    virtual scalar cl(scalar alfaRad, scalar reynolds, scalar mach)=0;
    virtual scalar cd(scalar alfaRad, scalar reynolds, scalar mach)=0;

    static Foam::autoPtr<Foam::airfoilModel> New(const word name, const dictionary& dict);

    //virtual scalar cm(scalar alfaRad, scalar reynolds, scalar mach);
protected:

    const word name_;

};

}

#endif /* PROPELLER_AIRFOIL_MODEL_AIRFOIL_MODEL */
