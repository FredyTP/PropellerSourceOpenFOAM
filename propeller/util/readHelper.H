#ifndef PROPELLER_UTIL_READ_HELPER
#define PROPELLER_UTIL_READ_HELPER
#include "dictionary.H"
#include "scalar.H"
#include "word.H"
#include "csvTable.H"
#include "Enum.H"
#include "autoPtr.H"
#include "List.H"
#include "label.H"

#include "InverseDistance.H"
#include "RegularInterpolation.H"
#include "LinearInterpolation.H"
#include "ClosestNeighbor.H"


namespace Foam
{

namespace util
{

autoPtr<RegularInterpolation<scalar,scalar,1>> NewInterpolationFromDict
(
    const dictionary& dict,
    word x_name,
    word y_name,
    bool convertToRad = false,
    bool enableCSV=false,
    const csvTable<scalar,word>* csv = nullptr
);

enum class interpolationMethod
{
    linear,
    closestNeighbor,
    inverseDistance
};

extern const Enum<interpolationMethod> interpolationNames;


template<class typeIn, class typeOu, label dim>
autoPtr<InterpolationTable<typeIn,typeOu,dim>> NewInterpolationFromRaw
(
    const dictionary& dict,
    const List<List<typeIn>>& xList,
    const List<typeOu>& yList
)
{
    typedef InterpolationTable<typeIn,typeOu,dim> baseType;
    typedef LinearInterpolation<typeIn,typeOu,dim> linearType;
    typedef ClosestNeighbor<typeIn,typeOu,dim> closestType;
    typedef InverseDistance<typeIn,typeOu,dim> inverseType;


    interpolationMethod method = interpolationMethod::linear;

    bool present = interpolationNames.readIfPresent("interpolation",dict,method);
    
    switch (method)
    {
    case interpolationMethod::linear:
    {
        autoPtr<baseType> ptr(autoPtr<linearType>::New());
        if(ptr->setRawData(xList,yList))
        {
            return ptr;
        }
        
        if(present)
        {
            FatalErrorInFunction
            <<"Data is not valid for linear Interpolation (structured)"
            <<exit(FatalError);
        }
    } 
    [[fallthrough]];
    case interpolationMethod::inverseDistance:
    {
        int nPoints = dict.getOrDefault("nPoints",-1);

        autoPtr<baseType> ptr(autoPtr<inverseType>::New());
        if(ptr->setRawData(xList,yList))
        {
            return ptr;
        }

        if(present)
        {
            FatalErrorInFunction
            <<"Data is not valid for inverseDistance Interpolation"
            <<exit(FatalError);
        }
    }
    [[fallthrough]];
    case interpolationMethod::closestNeighbor:
    {
        autoPtr<baseType> ptr(autoPtr<closestType>::New());
        if(ptr->setRawData(xList,yList))
        {
            return ptr;
        }
        FatalErrorInFunction
        <<"Data is not valid for closestType Interpolation"
        <<exit(FatalError);
        
    }
    [[fallthrough]];
    default:
        break;
    }

    return autoPtr<baseType>();
}

}
}

#endif /* PROPELLER_UTIL_READ_HELPER */
