#ifndef PROPELLER_ROTOR_MESH_ROTOR_GEOMETRY
#define PROPELLER_ROTOR_MESH_ROTOR_GEOMETRY
#include "scalar.H"
#include "vector.H"
#include "property.H"
#include "StringStream.H"
#define NO_RADIUS -1.0

namespace Foam
{
class rotorGeometry
{
private:
    util::property<scalar> innerRadius_;
    util::property<scalar> radius_;
    util::property<vector> center_;
    util::property<vector> direction_;
    util::property<vector> psiRef_;

public:

    rotorGeometry()
    {
        innerRadius_.addCondition(util::property<scalar>::isPositive);
        radius_.addCondition(util::property<scalar>::isGreaterThan(innerRadius_.get()));
        direction_.addCondition(isValidVector);
        psiRef_.addCondition(isValidVector);
    }
    

    util::property<scalar>& innerRadius() {return innerRadius_;}
    util::property<scalar>& radius() {return radius_;}
    util::property<vector>& center() {return center_;}
    util::property<vector>& direction() {return direction_;}
    util::property<vector>& psiRef() {return psiRef_;}


    const util::property<scalar>& innerRadius() const {return innerRadius_;}
    const util::property<scalar>& radius() const {return radius_;}
    const util::property<vector>& center() const {return center_;}
    const util::property<vector>& direction() const {return direction_;}
    const util::property<vector>& psiRef() const {return psiRef_;}
   
    /**
     * Read geometry presented data from dictionary
     * If a entry is preset, set the data and set the set flag to true
     * If a entry not present, dont modify the flag
    */
    bool readIfPresent(const dictionary& dict)
    {
        this->readIfPresent(dict,radius_,"radius");
        this->readIfPresent(dict,innerRadius_,"innerRadius");
        this->readIfPresent(dict,direction_,"direction");
        this->readIfPresent(dict,center_,"center");
        this->readIfPresent(dict,psiRef_,"psiRef");

        return isReady();
    }

    template<class pType>
    bool readIfPresent(const dictionary& dict,util::property<pType>& prop,word name)
    {
        pType temp;
        bool ok = dict.readIfPresent(name,temp);
        if(ok)
        {
            prop.set(temp);
        }
        return ok;
    }

    inline bool isReady()
    {
        return isSet() && isValid();
    }

    inline bool isValid()
    {
        return  innerRadius_.isValid() &&
                radius_.isValid()      &&
                direction_.isValid()   &&
                center_.isValid()      &&
                psiRef_.isValid();
    }

    inline bool isSet()
    {
        return  innerRadius_.isSet() && 
                radius_.isSet()      && 
                center_.isSet()      && 
                direction_.isSet()   && 
                psiRef_.isSet();
    }
    template<class pType>
    static string getFieldString(util::property<pType>& prop)
    {
        OStringStream field;
        /*if(!prop.isSet())
        {
            field << "unset";
        }
        else if(!prop.isValid())
        {
            field << "invalid("<<prop.get()<<")";
        }
        else
        {
            field<<prop.get();
        }*/
        return field.str();
    }
    friend Ostream& operator<<(Ostream& os, const rotorGeometry& r)
    {
        os.incrIndent();
        word innerRadius,radius,center,direction,psiRef;
        /*innerRadius <<"- innerRadius: "<<rotorGeometry::getFieldString(r.innerRadius())<<endl;
        indent(os) <<"- radius: "<<rotorGeometry::getFieldString(r.radius())<<endl;
        indent(os) <<"- center: "<<rotorGeometry::getFieldString(r.center())<<endl;
        indent(os) <<"- direction: "<<rotorGeometry::getFieldString(r.direction())<<endl;
        indent(os) <<"- psiRef: "<< rotorGeometry::getFieldString(r.psiRef())<<endl;*/

        os.decrIndent();
        return os;
    }
protected: 

    static bool isValidVector(const vector& v) { return magSqr(v) > VSMALL;}

    
};



}



#endif /* PROPELLER_ROTOR_MESH_ROTOR_GEOMETRY */
