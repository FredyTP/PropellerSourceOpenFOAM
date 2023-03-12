#ifndef PROPELLER_BLADE_MODEL_BLADE_SECTION
#define PROPELLER_BLADE_MODEL_BLADE_SECTION
#include "interpolated.H"
namespace Foam
{

class bladeSection
{
public:
    //- Airfoil of the section
    const airfoilModel* airfoil;

    //- Radius where the section is located in [m]
    scalar radius;

    //- Chort of section in [m]
    scalar chord;        

    //- Twisted angle of the section in [rad]
    scalar twistAngle;

    //- Sweep angle of the section in [rad]
    scalar sweepAngle;

    bladeSection(const airfoilModel* airfoil_,scalar radius_,scalar chord_,scalar twistAngle_, scalar sweepAngle_)
    : airfoil(airfoil_), radius(radius_),chord(chord_),twistAngle(twistAngle_),sweepAngle(sweepAngle_)
    {

    }  
    bladeSection()
    {

    }
};

class interpolatedBladeSection
: public interpolated<scalar,bladeSection>
{
public:

    interpolatedBladeSection()
    {

    }
    interpolatedBladeSection(interpolated<scalar,bladeSection>& other)
    {
        this->coeff = other.coefficients();
        this->nodes = other.points();
    }
    scalar radius()
    {
        return this->value([](scalar coeff, bladeSection b){return coeff*b.radius;});
    }
    scalar chord()
    {
        return this->value([](scalar coeff, bladeSection b){return coeff*b.chord;});
    }
    scalar twist()
    {
        return this->value([](scalar coeff, bladeSection b){return coeff*b.twistAngle;});
    }
    scalar sweep()
    {
        return this->value([](scalar coeff, bladeSection b){return coeff*b.sweepAngle;});
    }
    scalar cl(scalar alpha, scalar reynolds, scalar mach)
    {
        return this->value([=](scalar coeff, bladeSection b){return coeff*b.airfoil->cl(alpha,reynolds,mach);});
    }
    scalar cd(scalar alpha, scalar reynolds, scalar mach)
    {
        return this->value([=](scalar coeff, bladeSection b){return coeff*b.airfoil->cd(alpha,reynolds,mach);});
    }
};

}

#endif /* PROPELLER_BLADE_MODEL_BLADE_SECTION */
