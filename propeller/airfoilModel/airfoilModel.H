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
                word name,
                const dictionary& dict
            ),
            (name, dict)
        );

    airfoilModel();
    airfoilModel(word name);
    airfoilModel(word name, const dictionary& dict);
    
    virtual ~airfoilModel() = default;
    
    const word& airfoilName() const {return name_;}
    
    virtual scalar cl(scalar alfaRad, scalar reynolds, scalar mach) const = 0;
    virtual scalar cd(scalar alfaRad, scalar reynolds, scalar mach) const = 0;

    static Foam::autoPtr<Foam::airfoilModel> New(const word name, const dictionary& dict);

    void exportAirfoil();
    void writeAirfoil(fileName path, scalar alfaBegin, scalar alfaEnd, label nAlfa, const List<scalar>& Reyn, const List<scalar>& Mach);

    //virtual scalar cm(scalar alfaRad, scalar reynolds, scalar mach);
private:

    word name_;

    //Export properties
    bool export_;
    fileName path_;
    scalar aoaBegin_;
    scalar aoaEnd_;
    label nAoa_;
    List<scalar> Re_;
    List<scalar> Ma_;

};



}

#endif /* PROPELLER_AIRFOIL_MODEL_AIRFOIL_MODEL */
