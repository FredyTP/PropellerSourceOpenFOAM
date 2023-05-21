#ifndef PROPELLER_VELOCITY_SAMPLER_TEST_TEMPLATED
#define PROPELLER_VELOCITY_SAMPLER_TEST_TEMPLATED

#include "dictionary.H"
#include "autoPtr.H"
#include "rotorGrid.H"
#include "rotorFvMeshSel.H"
#include "runTimeSelectionTables.H"

#include "addToRunTimeSelectionTable.H"


namespace Foam
{
template<class T>
class fieldSampler
{
public:
    
    TypeName("fieldSampler");

    fieldSampler(const dictionary& dict) {}
    declareRunTimeSelectionTable
    (
        autoPtr,
        fieldSampler,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );
    static fieldSampler New(const dictionary& dict);
     
    virtual ~fieldSampler() = default;
};

typedef fieldSampler<scalar> scalarFieldSampler;


defineTemplateTypeNameWithName(fieldSampler<scalar>,"scalarFieldSampler");
defineTemplatedRunTimeSelectionTable(fieldSampler, dictionary,scalar);


template<class fieldType>
class domainSampler : public fieldSampler<fieldType>
{
public:
    TypeName("domainSampler");
    domainSampler(const dictionary& dict) : fieldSampler<fieldType>(dict) {Info<<"CREATED DOMAIN SAMPLER OMG"<<endl;}

};

defineTemplateTypeNameWithName(domainSampler<scalar>,"scalarDomainSampler");
addTemplatedToRunTimeSelectionTable(fieldSampler,domainSampler,scalar,dictionary);

template <class T>
inline fieldSampler<T> fieldSampler<T>::New(const dictionary &dict)
{
    //Get model Type name (Ex: fixedVelocity) 
        //From type key from dictionary (propellerModel)
        const word modelType("scalarDomainSampler"); 
        Info<<endl;
        Info<< "Selecting " << typeName << " " << modelType << endl;

        //Find class contructor in tables
        auto* ctorPtr = dictionaryConstructorTable(modelType);

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                dict,
                typeName,
                modelType,
                *dictionaryConstructorTablePtr_
            ) << exit(FatalIOError);
        }

        return autoPtr<fieldSampler<T>>(ctorPtr(dict));

}

}

#endif /* PROPELLER_VELOCITY_SAMPLER_TEST_TEMPLATED */
