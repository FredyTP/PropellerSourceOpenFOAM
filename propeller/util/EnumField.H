#ifndef PROPELLER_UTIL_ENUM_FIELD
#define PROPELLER_UTIL_ENUM_FIELD
#include "scalarField.H"
#include "label.H"
#include "Enum.H"

namespace Foam
{
namespace util
{


template<class enumType,label N>
class EnumField : public scalarField
{
public:
    EnumField()
    : scalarField(N,0)
    {

    }
    EnumField(scalar val)
        : scalarField(N,val)
    {

    }
    List<enumType> readIfPresent(const dictionary& dict, const Enum<enumType>& enumNames)
    {
        const List<word>& names = enumNames.names();
        const List<int>& values = enumNames.values();
        List<enumType> found;

        forAll(names,i)
        {
            if(dict.readIfPresent(names[i], scalarField::operator[](values[i])))
            {
                Info<<"i: "<<i<<endl;
                Info<<"name: "<<names[i]<<endl;
                Info<<"enum: "<<values[i]<<endl;
                found.append(static_cast<enumType>(values[i]));
            }
        }
        return found;
    }
    void readOrDefault(const dictionary& dict, const Enum<enumType>& enumNames, scalar defaultVal)
    {
        const List<word>& names = enumNames.names();
        const List<int>& values = enumNames.values();

        forAll(names,i)
        {
            scalarField::operator[](values[i]) =  dict.getOrDefault<scalar>(names[i], defaultVal);
        }
    }
    void set(const scalarField& values, const List<enumType>& pos)
    { 
        if(values.size()!= pos.size())
        {
            return;
        }
        
        forAll(values,i)
        {
            (*this)[pos[i]]=values[i];
        }
    }

    void get(scalarField& values, const List<enumType>& pos) const
    {
        if(values.size()!= pos.size())
        {
            values.resize(pos.size());
        }
        
        forAll(values,i)
        {
            values[i] = (*this)[pos[i]];
        }
    }
    scalar& operator[](enumType val)
    {
        return scalarField::operator[](static_cast<label>(val));
    }
    scalar operator[](enumType val) const
    {
        return scalarField::operator[](static_cast<label>(val));
    }
};



}
}

#endif /* PROPELLER_UTIL_ENUM_FIELD */
