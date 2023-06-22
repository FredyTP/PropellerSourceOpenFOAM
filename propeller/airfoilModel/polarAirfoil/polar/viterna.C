#include "viterna.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
namespace Foam
{
    //set and define the type name "typeName"
    defineTypeNameAndDebug(viterna,0);

    //Add to run time table to dynamically select the class based on
    //the dictionary airfoilModel atribute

    addToRunTimeSelectionTable(polar,viterna,dictionary);


viterna::viterna(bool cubicSpline, List<scalar> alpha, List<scalar> cl, List<scalar> cd, scalar Re, scalar Ma, bool isRadian)
 : polar(cubicSpline,alpha,cl,cd,Re,Ma,isRadian)
{
    setParameters();
}

viterna::viterna(bool cubicSpline, const fileName filename, scalar Re, scalar Ma, bool isRadian)
: polar(cubicSpline,filename,Re,Ma,isRadian)
{
    setParameters();
}

scalar viterna::cl(scalar alpha) const
{
    using namespace Foam::constant::mathematical;
    //alpha in rads between [-pi, pi]
    if(alpha>=data_.alpha_stall_neg && alpha<=data_.alpha_stall_pos)
    {
        return polar::cl(alpha);
    }
    else if(alpha>data_.alpha_stall_pos && alpha <=piByTwo) // [stall, 90º]
    {
        return viterna_cl(A1_p,A2_p,alpha);
    }
    else if((alpha>piByTwo && alpha<=pi) || (alpha>=-pi && alpha<-piByTwo)) // [90 , 180] o [-180, -90]
    {
        return flat_plate_cl(data_.cd_max,alpha);
    }
    else // [-90, negative stall]
    {
        return -viterna_cl(A1_n,A2_n,-alpha); //antisimetric transformation -alfa, -cl
    }
        
}

scalar viterna::cd(scalar alpha) const
{
    using namespace Foam::constant::mathematical;
    //alpha in rads between [-pi, pi]
    if(alpha>=data_.alpha_stall_neg && alpha<=data_.alpha_stall_pos)
    {
        return polar::cd(alpha);
    }
    else if(alpha>data_.alpha_stall_pos && alpha <=piByTwo) // [stall, 90º]
    {
        return viterna_cd(B1_p,B2_p,alpha);
    }
    else if((alpha>piByTwo && alpha<=pi) || (alpha>=-pi && alpha<-piByTwo)) // [90 , 180] o [-180, -90]
    {
        return flat_plate_cd(data_.cd_max,alpha);
    }
    else // [-90, negative stall]
    {
        return viterna_cd(B1_n,B2_n,-alpha); //antisimetric transformation -alfa, -cl
    }
}

scalar viterna::calcA1(scalar cd_max) const
{
    return cd_max/2;
}
scalar viterna::calcA2(scalar cd_max, scalar cl_stall, scalar alpha_stall) const
{
    return (cl_stall - cd_max*sin(alpha_stall)*cos(alpha_stall))
            *sin(alpha_stall)/pow(cos(alpha_stall),2);
}
scalar viterna::calcB1(scalar cd_max) const
{
    return cd_max;
}
scalar viterna::calcB2(scalar cd_max, scalar cd_stall, scalar alpha_stall) const
{
    return (cd_stall - cd_max * pow(sin(alpha_stall),2)) / (cos(alpha_stall));
}

scalar viterna::viterna_cl(scalar A1, scalar A2, scalar alpha) const
{
    return A1*sin(2*alpha)+A2*(pow(cos(alpha),2))/sin(alpha);
}

scalar viterna::viterna_cd(scalar B1, scalar B2, scalar alpha) const
{
    return B1*pow(sin(alpha),2) + B2*cos(alpha);
}

scalar viterna::flat_plate_cl(scalar cl_max_flat_plate, scalar alpha) const
{
    return 2*std::abs(cl_max_flat_plate)*sin(alpha)*cos(alpha);
}
scalar viterna::flat_plate_cd(scalar cd_max_flat_plate, scalar alpha) const
{
    return cd_max_flat_plate * pow(sin(alpha),2);
}
void viterna::setParameters()
{

    //Set data from polar
    scalar AR=10; //Usually works (?)
    data_.alpha_stall_pos = alpha_max; //rad
    data_.alpha_stall_neg = alpha_min; //rad 
    data_.cl_stall_pos = cl_alpha_max;
    data_. cl_stall_neg = cl_alpha_min;
    data_.cd_stall_pos=cd_alpha_max;
    data_.cd_stall_neg = cd_alpha_min;
    data_.cd_max = 1.11 + 0.018*AR;

    //Calculate extrapolation parameters
    A1_p = calcA1(data_.cd_max);
    A1_n = calcA1(data_.cd_max);

    A2_p = calcA2(data_.cd_max,std::abs(data_.cl_stall_pos),std::abs(data_.alpha_stall_pos));
    A2_n = calcA2(data_.cd_max,std::abs(data_.cl_stall_neg),std::abs(data_.alpha_stall_neg));

    B1_p = calcB1(data_.cd_max);
    B1_n = calcB1(data_.cd_max);

    B2_p = calcB2(data_.cd_max,std::abs(data_.cd_stall_pos),std::abs(data_.alpha_stall_pos));
    B2_n = calcB2(data_.cd_max,std::abs(data_.cd_stall_neg),std::abs(data_.alpha_stall_neg));
}

bool viterna::valid()
{
    bool valid = true;

    if(alpha_min >= 0 || alpha_max <= 0)
    {
        valid = false;
    }

    return (valid & polar::valid());
}

}

