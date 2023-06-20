#ifndef PROPELLER_UTIL_OUTPUT_FIELDS
#define PROPELLER_UTIL_OUTPUT_FIELDS

#include "volFields.H"
#include "Enum.H"
#include "HashTable.H"

namespace Foam
{
namespace util
{
    


/**
 * @brief Class that holds volFields to be used as output,
 * being able to enable or disable them
 * 
 * @tparam fieldType 
 * @tparam enumType 
 * @tparam N 
 */
template<class fieldType,class enumType, label N>
class outputFields
{
public:
    outputFields
    (
        word ownerName,
        const dictionary& dict,
        const fvMesh& mesh,
        std::initializer_list<std::pair<enumType, const char*>> enumNames,
        std::initializer_list<std::pair<label, dimensionSet>> enumDimensions
    );

    void read(const dictionary& dict);
    void initFields();
    void resetFields();

    /**
     * @brief Set the field value to specified cell 
     * 
     * @param field 
     * @param value 
     * @param cell 
     */
    void setValue(enumType field, const fieldType& value, label cell);

    /**
     * @brief Set the field value to specified cells 
     * 
     * @param field 
     * @param value 
     * @param cells 
     */
    void setValue(enumType field, const fieldType& value, const List<label>& cells);

    /**
     * @brief Set the field value to specified cells with weighted values
     * 
     * @param field 
     * @param value 
     * @param cells 
     * @param weights 
     */
    void setValue(enumType field, const fieldType& value, const List<label>& cells, const List<scalar>& weights);


    /**
     * @brief Apply values to enabled fields
     * 
     * @param values 
     * @param cells 
     */
    void setValues(const FixedList<fieldType,N>& values, const List<label>& cells);


    
private:
    word ownerName_;
    const fvMesh& mesh_;
    Enum<enumType> names_;
    HashTable<dimensionSet,label> dimensions_;
    List<bool> active_;

    PtrList<VolumeField<fieldType>> fields_;
};

template <class fieldType, class enumType, label N>
inline outputFields<fieldType, enumType, N>::outputFields
(   
    word ownerName,
    const dictionary& dict,
    const fvMesh& mesh,
    std::initializer_list<std::pair<enumType, const char*>> enumNames,
    std::initializer_list<std::pair<label, dimensionSet>> enumDimensions
) : 
    ownerName_(ownerName),
    mesh_(mesh),
    names_(enumNames),
    dimensions_(enumDimensions),
    active_(N,false),
    fields_(N)
{
    this->read(dict);
    this->initFields();
}

template <class fieldType, class enumType, label N>
inline void outputFields<fieldType, enumType, N>::read(const dictionary &dict)
{
    List<word> enabledNames = dict.getOrDefault<List<word>>("outputFields", List<word>());
    forAll(enabledNames,i)
    {
        const word& name = enabledNames[i];
        label idx = names_.find(name);
        if( idx != -1)
        {
            active_[idx] = true;
        }
    }
}

template <class fieldType, class enumType, label N>
inline void outputFields<fieldType, enumType, N>::initFields()
{
    forAll(active_,i)
    {
        if(!active_[i])
        {
            continue;
        }

        VolumeField<fieldType>* fieldPtr = new VolumeField<fieldType>
        (
            IOobject
            (
                ownerName_ +":" +names_.get(static_cast<enumType>(i)),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensioned<fieldType>(dimless, Zero)
        );
        fields_.set(i,fieldPtr);
    }
}

template <class fieldType, class enumType, label N>
inline void outputFields<fieldType, enumType, N>::resetFields()
{
    forAll(active_,i)
    {
        if(!active_[i])
        {
            continue;
        }
        fields_[i].ref(false) = dimensioned<fieldType>(dimensions_.at(i), Zero);
    }
}

template <class fieldType, class enumType, label N>
inline void outputFields<fieldType, enumType, N>::setValue(enumType field, const fieldType &value, label cell)
{
    label idx = static_cast<label>(field);
    if(active_[idx])
    {
        fields_[idx].ref(false)[cell]=value;
    }
}

template <class fieldType, class enumType, label N>
inline void outputFields<fieldType, enumType, N>::setValue(enumType field, const fieldType &value, const List<label> &cells)
{
    label idx = static_cast<label>(field);
    if(active_[idx])
    {
        forAll(cells,i)
        {
            fields_[idx].ref(false)[cells[i]]=value;
        }
    }
}

template <class fieldType, class enumType, label N>
inline void outputFields<fieldType, enumType, N>::setValue(enumType field, const fieldType &value, const List<label> &cells, const List<scalar> &weights)
{
    label idx = static_cast<label>(field);
    if(active_[idx])
    {
        forAll(cells,i)
        {
            fields_[idx].ref(false)[cells[i]]=weights[i]*value;
        }
    }
}
template <class fieldType, class enumType, label N>
inline void outputFields<fieldType, enumType, N>::setValues(const FixedList<fieldType,N> &values, const List<label> &cells)
{
    forAll(values,i)
    {
        if(active_[i])
        {
            forAll(cells,i)
            {
                fields_[i].ref(false)[cells[i]]=values[i];
            }
        }
    }
}
}
}


#endif /* PROPELLER_UTIL_OUTPUT_FIELDS */


