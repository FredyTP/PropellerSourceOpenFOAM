#include "viterna.H"

namespace Foam
{

viterna::viterna(const viterna_data &data)
{
    data_ = viterna_data;
    A1_p = calcA1(data_.cd_max);
    A1_n = calcA1(data_.cd_max);

    A2_p = calcA2(data_.cd_max,abs(data_.cl_stall_pos),abs(data_.alpha_stall_pos));
    A2_n = calcA2(data_.cd_max,abs(data_.cl_stall_neg),abs(data_.alpha_stall_neg));

    B1_p = calcB1(data_.cd_max);
    B1_n = calcB1(data_.cd_max);

    B2_p = calcB2(data_.cd_max,abs(data_.cd_stall_pos),abs(data_.alpha_stall_pos));
    B2_n = calcB2(data_.cd_max,abs(data_.cd_stall_neg),abs(data_.alpha_stall_neg));
}

scalar viterna::cl(scalar alpha)
{
    return scalar();
}

scalar viterna::cd(scalar alpha)
{
    return scalar();
}

scalar viterna::calcA1(scalar cd_max)
{
    return cd_max/2;
}
scalar viterna::calcA2(scalar cd_max, scalar cl_stall, scalar alpha_stall)
{
    return (cl_stall - cd_max*sin(alpha_stall)*cos(alpha_stall))
            *sin(alpha_stall)/(cos(alpha_stall)^2);
}
scalar viterna::calcB1(scalar cd_max)
{
    return cd_max;
}
scalar viterna::calcB2(scalar cd_max, scalar cd_stall, scalar alpha_stall){
    return (cd_stall - cd_max * sin(alpha_stall) * *2) / (cos(alpha_stall))}

scalar viterna::viterna_cl(scalar A1, scalar A2, scalar alpha)
{
    return A1*sin(2*alpha)+A2*(cos(alpha)^2)/sin(alpha);
}

scalar viterna::viterna_cd(scalar B1, scalar B2, scalar alpha)
{
    return B1*sin(alpha)^2 + B2*cos(alpha);
}

scalar viterna::flat_plate_cl(scalar cl_max_flat_plate, scalar alpha)
{
    return 2*abs(cl_max_flat_plate)*sin(alpha)*cos(alpha);
}
scalar viterna::flat_plate_cd(scalar cd_max_flat_plate, scalar alpha)
{
    return 2*cd_max_flat_plate * sin(alpha)^2;
}
}