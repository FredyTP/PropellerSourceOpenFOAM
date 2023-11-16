#ifndef PROPELLER_AIRFOIL_MODEL_AIRFOIL_MODEL
#define PROPELLER_AIRFOIL_MODEL_AIRFOIL_MODEL

#include "dictionary.H"
#include "scalar.H"
#include "autoPtr.H"
#include "Interpolated.H"
namespace Foam
{


/**
 * @brief This abstract class is the base for all derived airfoils.
 * Contains virtual function for cl and cd and automatization of the export.
 * Contains runtimetable.
 * 
 */
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
    
    /**
     * @brief Virtual method for cl
     * 
     * @param alfaRad 
     * @param reynolds 
     * @param mach 
     * @return scalar 
     */
    virtual scalar cl(scalar alfaRad, scalar reynolds, scalar mach) const = 0;
    /**
     * @brief Virtual method for cd
     * 
     * @param alfaRad 
     * @param reynolds 
     * @param mach 
     * @return scalar 
     */
    virtual scalar cd(scalar alfaRad, scalar reynolds, scalar mach) const = 0;

    /**
     * @brief Creates new airfoil instance from runtime table and dictionary
     * 
     * @param name 
     * @param dict 
     * @return Foam::autoPtr<Foam::airfoilModel> 
     */
    static Foam::autoPtr<Foam::airfoilModel> New(const word name, const dictionary& dict);

    /**
     * @brief Exports current airfoil with configured params, enabling split or joint files
     * 
     */
    void exportAirfoil();

    /**
     * @brief Write airfoil data for configured params in one file
     * 
     * @param path 
     * @param alfaBegin 
     * @param alfaEnd 
     * @param nAlfa 
     * @param Reyn 
     * @param Mach 
     */
    void writeAirfoil(fileName path, scalar alfaBegin, scalar alfaEnd, label nAlfa, const List<scalar>& Reyn, const List<scalar>& Mach);

    //virtual scalar cm(scalar alfaRad, scalar reynolds, scalar mach);
private:

    //- Airfoil name, for id
    word name_;

    //- Export properties
    bool export_;
    bool join_;
    fileName path_;
    scalar aoaBegin_;
    scalar aoaEnd_;
    label nAoa_;
    List<scalar> Re_;
    List<scalar> Ma_;

};



}

#endif /* PROPELLER_AIRFOIL_MODEL_AIRFOIL_MODEL */
