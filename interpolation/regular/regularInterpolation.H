#ifndef PROPELLER_INTERPOLATION_REGULAR_REGULAR_INTERPOLATION
#define PROPELLER_INTERPOLATION_REGULAR_REGULAR_INTERPOLATION

#include "List.H"
#include "scalar.H"
#include "interpolationTable.H"
namespace Foam
{

/**
 * It is possible to make a N-dimesional regular interpolation
 * using a 1D array with a N-dimesional index mapping
 * index = N*i + j;
 * index = Nx*Ny*i + Nx*j + k; and so on...
 * 
 * the number of coefficients can be "easily" obtained in matrix form
*/
template<label dimension>
class regularInterpolation 
: public interpolationTable<dimension>
{
protected:
public:

    virtual scalar interpolate(FixedList<scalar,dimension> input)=0;

};

class regularInterpolation1D : public interpolationTable<1>
{

};

class regularInterpolation2D : public interpolationTable<2>
{

};

class regularInterpolation3D : public interpolationTable<3>
{

};


}



#include "regularInterpolation.C"

#endif /* PROPELLER_INTERPOLATION_REGULAR_REGULAR_INTERPOLATION */
