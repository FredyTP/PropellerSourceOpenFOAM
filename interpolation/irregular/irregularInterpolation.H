#ifndef PROPELLER_INTERPOLATION_IRREGULAR_IRREGULAR_INTERPOLATION
#define PROPELLER_INTERPOLATION_IRREGULAR_IRREGULAR_INTERPOLATION

#include "List.H"
#include "scalar.H"
namespace Foam
{

template<unsigned int dimension>
class irregularInterpolation
{
protected:
    const List<FixedList<scalar,dimension>>& input_;
    const List<scalar>& output_;
public:

    irregularInterpolation
    (
        const List<FixedList<scalar,dimension>>& input,
        const List<scalar>& output
    );

    

    virtual scalar interpolate(FixedList<scalar,dimension> input)=0;

    //Specialized operators () for variable size input
    scalar operator()(scalar x0);
    scalar operator()(scalar x0, scalar x1);
    scalar operator()(scalar x0, scalar x1,scalar x2);

};

template<>
scalar irregularInterpolation<1>::operator()(scalar x0)
{
    return this->interpolate({x0});
}
template<>
scalar irregularInterpolation<2>::operator()(scalar x0, scalar x1)
{
    return this->interpolate({x0,x1});
}
template<>
scalar irregularInterpolation<3>::operator()(scalar x0, scalar x1,scalar x2)
{
    return this->interpolate({x0,x1,x2});
}


}

#include "irregularInterpolation.C"


#endif /* PROPELLER_INTERPOLATION_IRREGULAR_IRREGULAR_INTERPOLATION */
