#include "interpolatedAirfoil.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    //set and define the type name "typeName"
    defineTypeNameAndDebug(interpolatedAirfoil,0);

    //Add to run time table to dynamically select the class based on
    //the dictionary airfoilModel atribute

    addToRunTimeSelectionTable(airfoilModel,interpolatedAirfoil,dictionary);

interpolatedAirfoil::interpolatedAirfoil()
{
}

interpolatedAirfoil::interpolatedAirfoil(word name,const dictionary& dict)
    : airfoilModel(name)
{
    List<scalar> coefficients = dict.get<List<scalar>>("coefficients");
    airfoilNames_ = dict.get<List<word>>("airfoils");

    if(coefficients.size()!= airfoilNames_.size())
    {
        FatalErrorInFunction<<"Coefficients size is diferent than airfoils: "
        <<coefficients.size()<<" != "<< airfoilNames_.size()
        <<exit(FatalError);
    }

    this->coeff = coefficients;
}

interpolatedAirfoil::interpolatedAirfoil(interpolated<scalar, const airfoilModel *> &other)
{
    this->coeff = other.coefficients();
    this->nodes = other.points();
}

void interpolatedAirfoil::build(const airfoilModelList& airfoilList)
{
    this->nodes.resize(airfoilNames_.size());
    forAll(airfoilNames_,i)
    {
        this->nodes[i] = airfoilList.getAirfoil(airfoilNames_[i]);
        if(this->nodes[i]==nullptr)
        {
            FatalErrorInFunction
            <<"Cannot find airfoil: "<<airfoilNames_[i]
            <<exit(FatalError);
        }

    }
}

scalar Foam::interpolatedAirfoil::cl(scalar alfaRad, scalar reynolds, scalar mach) const
{
    return this->value<scalar>([=](scalar coeff, const airfoilModel* air){return coeff*air->cl(alfaRad,reynolds,mach);});
}

scalar interpolatedAirfoil::cd(scalar alfaRad, scalar reynolds, scalar mach) const
{
    return this->value<scalar>([=](scalar coeff, const airfoilModel* air){return coeff*air->cd(alfaRad,reynolds,mach);});
}

}

