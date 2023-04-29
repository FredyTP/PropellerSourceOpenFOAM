#include "interpolatedAirfoil.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    //set and define the type name "typeName"
    //defineTypeNameAndDebug(interpolatedAirfoil,0);

    //Add to run time table to dynamically select the class based on
    //the dictionary airfoilModel atribute

    //addToRunTimeSelectionTable(airfoilModel,interpolatedAirfoil,dictionary);

interpolatedAirfoil::interpolatedAirfoil()
{
}

interpolatedAirfoil::interpolatedAirfoil(interpolated<scalar, const airfoilModel *> &other)
{
    this->coeff = other.coefficients();
    this->nodes = other.points();
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