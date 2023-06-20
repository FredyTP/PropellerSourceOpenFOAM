#ifndef PROPELLER_UTIL_ENUM_FIELD
#define PROPELLER_UTIL_ENUM_FIELD
#include "scalarField.H"
#include "label.H"
#include "Enum.H"

namespace Foam
{
namespace util
{

/**
 * @brief Class that stores, manipulate and give access to an array of values of a Enum field
 * 
 * @tparam enumType 
 * @tparam N 
 */
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

    /**
     * @brief Read from dictionary the found enums and the value specified and returns the found enums
     * 
     * @param dict 
     * @param enumNames enumName to enum map
     * @return List<enumType> found enums in found order
     */
    List<enumType> readIfPresent(const dictionary& dict, const Enum<enumType>& enumNames)
    {
        const List<word>& names = enumNames.names();
        const List<int>& values = enumNames.values();
        List<enumType> found;

        forAll(names,i)
        {
            if(dict.readIfPresent(names[i], scalarField::operator[](values[i])))
            {
                found.append(static_cast<enumType>(values[i]));
            }
        }
        return found;
    }
    /**
     * @brief Read the enum values from dict. Assign default value if entry not found.
     * 
     * @param dict 
     * @param enumNames enumName to enum map
     * @param defaultVal default value assigned to missing entries
     */
    void readOrDefault(const dictionary& dict, const Enum<enumType>& enumNames, scalar defaultVal)
    {
        const List<word>& names = enumNames.names();
        const List<int>& values = enumNames.values();

        forAll(names,i)
        {
            scalarField::operator[](values[i]) =  dict.getOrDefault<scalar>(names[i], defaultVal);
        }
    }

    /**
     * @brief Set the internal field values from the input values and enum index
     * 
     * @param values setted values
     * @param pos enum index of the values
     */
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

    /**
     * @brief Get in the values list the selected enum values in specified order
     * 
     * @param values list to get filled
     * @param pos selected enum values to obtain
     */
    void get(scalarField& values, const List<enumType>& pos) const
    {
        values.resize(pos.size());
        
        forAll(values,i)
        {
            values[i] = (*this)[pos[i]];
        }
    }

    /**
     * @brief Overloaded [] operator to acces or modify the internal field throw the enum
     * 
     * @param val 
     * @return scalar& 
     */
    scalar& operator[](enumType val)
    {
        return scalarField::operator[](static_cast<label>(val));
    }
    /**
     * @brief Overloaded const [] operator to acces (only) the internal field throw the enum
     * 
     * @param val 
     * @return scalar 
     */
    scalar operator[](enumType val) const
    {
        return scalarField::operator[](static_cast<label>(val));
    }
};



}
}

#endif /* PROPELLER_UTIL_ENUM_FIELD */
