#include "fixedValue.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam
{
 

defineTemplateTypeNameWithName(fixedValue<vector>,"fixedValue");
addTemplatedToRunTimeSelectionTable(diskSampler,fixedValue,vector,dictionary);

defineTemplateTypeNameWithName(fixedValue<scalar>,"fixedValue");
addTemplatedToRunTimeSelectionTable(diskSampler,fixedValue,scalar,dictionary);



template<class fType>
fixedValue<fType>::fixedValue(const dictionary& dict,const rotorGrid* rGrid,const rotorFvMeshSel* rMesh)
    : diskSampler<fType>(rGrid,rMesh)
{
    this->read(dict);
}

template<class fType>
bool fixedValue<fType>::read(const dictionary &dict)
{
    Info.stream().incrIndent();

    bool normal = dict.getOrDefault<bool>("normal","false");
    bool ok=true;
    
    if(normal)
    {
        scalar speed;
        ok &=dict.readEntry("value",speed);
        //Positive speed towards the disk
        //fieldValue_ = - speed * this->rGrid_->geometry().direction().get();
        indent(Info)<< "- Normal to rotor speed: "<<speed<<endl;
    }
    else
    {
        ok &=dict.readEntry("value",fieldValue_);
    }

    indent(Info)<< "- Rotor input velocity: "<<fieldValue_<<endl;

    forAll(this->sampledField_,i)
    {
        this->sampledField_[i]=fieldValue_;
    }

    Info.stream().decrIndent();
    return ok;
    
}
template<class fType>
const Field<fType>& fixedValue<fType>::sampleField(const GeometricField<fType, fvPatchField, volMesh>& U)
{
    return this->sampledField_;
}

}
