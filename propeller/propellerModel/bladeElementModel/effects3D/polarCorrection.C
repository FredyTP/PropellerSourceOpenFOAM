#include "polarCorrection.H"

namespace Foam
{
namespace effects3D
{
const Enum
<
    polarCorrection::model
>
polarCorrection::modelNames_
({
    {model::None, "none"},
    {model::Snel, "Snel"},
    {model::SnelPumping, "SnelPumping"}
});

polarCorrection::polarCorrection(const dictionary &dict)
{
    usedModel_ = modelNames_.getOrDefault("polarCorrection",dict,model::None);
}

void polarCorrection::correct(scalar &cl3d, scalar cl2d, scalar alpha, scalar chord, scalar radius, scalar maxRadius, scalar angularVelocity, scalar Veff) const
{
    switch (usedModel_)
    {
    case model::None :
        cl3d=cl2d;
        break;
    case model::Snel :
        polarCorrection::Snel(cl3d,cl2d,alpha,chord,radius,maxRadius);
        break;
    case model::SnelPumping :
        polarCorrection::SnelPumping(cl3d,cl2d,alpha,chord,radius,maxRadius,angularVelocity,Veff);
        break;
    default:
        break;
    }
}

void polarCorrection::Snel(scalar &cl3d, scalar cl2d, scalar alpha, scalar chord, scalar radius, scalar maxRadius)
{
    //Only correct for alpha between 0 and 50 deg
    if(alpha < 0.0 || alpha > deg50inRad || cl2d <= VSMALL)
    {
        cl3d=cl2d;
        return;
    }
    //Only correct for radius < 80%
    if(radius/maxRadius > 0.8)
    {
        cl3d=cl2d;
        return;
    }

    scalar cl_pot = constant::mathematical::twoPi*alpha;
    if(cl_pot<cl2d)
    {
        cl3d=cl2d;
    }

    scalar Dcl = 3.1*pow(chord/radius,2)*(cl_pot-cl2d);
    scalar highAlphaFactor = 1.0;
    if(alpha>deg30inRad) 
    {
        highAlphaFactor = 1.0 - (alpha-deg30inRad)/(deg50inRad-deg30inRad);
    } 
    
    cl3d = cl2d + highAlphaFactor * Dcl;
}

void polarCorrection::SnelPumping(scalar &cl3d, scalar cl2d, scalar alpha, scalar chord, scalar radius, scalar maxRadius, scalar angularVelocity, scalar Veff)
{
    //Only correct for alpha between 0 and 50 deg
    if(alpha < 0.0 || alpha >deg50inRad || radius < VSMALL || Veff < VSMALL)
    {
        cl3d=cl2d;
        return;
    }
    //Only correct for radius < 80%
    if(radius/maxRadius>0.8)
    {
        cl3d=cl2d;
        return;
    }
    scalar cl_pot = constant::mathematical::twoPi*alpha;
    if(cl_pot<cl2d)
    {
        cl3d=cl2d;
    }

    scalar Dcl = 3.1*pow((angularVelocity*radius*chord)/(Veff*radius),2)*(cl_pot-cl2d);
    scalar highAlphaFactor = 1.0;
    if(alpha>deg30inRad) 
    {
        highAlphaFactor = 1.0 - (alpha-deg30inRad)/(deg20inRad);
    } 
    
    cl3d = cl2d + highAlphaFactor * Dcl;
}
}
}
