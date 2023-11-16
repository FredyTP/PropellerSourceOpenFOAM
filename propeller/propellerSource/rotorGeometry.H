#ifndef PROPELLER_PROPELLER_SOURCE_ROTOR_GEOMETRY
#define PROPELLER_PROPELLER_SOURCE_ROTOR_GEOMETRY
#include "scalar.H"
#include "vector.H"
#include "property.H"
#include "cartesianCS.H"
#include "cylindricalCS.H"

namespace Foam
{
/**
 * @brief Contains the definition of the rotor geometry and coordinate systems
 * 
 */
class rotorGeometry
{
private:
    util::property<scalar> innerRadius_;
    util::property<scalar> radius_;
    util::property<vector> center_;
    util::property<vector> direction_;
    util::property<vector> psiRef_;

    coordSystem::cartesian cartesianCS_;
    coordSystem::cylindrical cylindricalCS_;
    coordSystem::cylindrical cartesianToCylindrical_;


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
    
    friend Ostream& operator<<(Ostream& os, const rotorGeometry& r)
    {
        os.incrIndent();
        indent(os) <<"- innerRadius: "<<r.innerRadius()<<endl;
        indent(os) <<"- radius: "<<r.radius()<<endl;
        indent(os) <<"- center: "<<r.center()<<endl;
        indent(os) <<"- direction: "<<r.direction()<<endl;
        indent(os) <<"- psiRef: "<< r.psiRef()<<endl;

        os.decrIndent();
        return os;
    }

    /**
     * @brief Builds the 3 different coordinate systems required
     * to transform from global to local points
     * 
     */
    void buildCS()
    {
        cylindricalCS_ =  coordSystem::cylindrical // local Cartesian to cylindrical
        (
            center_.get(),    // centerd to local
            direction_.get(), // z-axis
            psiRef_.get()     // x-axis
        );


        cartesianCS_ = coordSystem::cartesian(
            center_.get(),    // centerd to local
            direction_.get(), // z-axis
            psiRef_.get()     // x-axis
        );

        cartesianToCylindrical_ = coordSystem::cylindrical(
            vector(0, 0, 0), // same center
            vector(0, 0, 1), // same z-axis
            vector(1, 0, 0)  // same x-axis
        );
    }

    /**
     * @brief Return cylindrical Coordinate system
     * to convert from global points (x y z) to local cilindrical points (r theta z)
     * 
     * @return coordSystem::cylindrical 
     */
    const coordSystem::cylindrical& cylindricalCS() const
    {
        return cylindricalCS_;
    }
    /**
     * @brief Return cartesian Coordinate system
     * to convert from global points (x y z) to local cartesian points (x y z)
     * 
     * @return coordSystem::cylindrical 
     */
    const coordSystem::cartesian& cartesianCS() const 
    {
        return cartesianCS_;
    }
    /**
     * @brief Return cylindrical Coordinate system
     * to convert from local cartesian points (x y z) to local cilindrical points (r- theta - z)
     * 
     * @return coordSystem::cylindrical 
     */
    const coordSystem::cylindrical& cartesianToCylindrical() const
    {
        return cartesianToCylindrical_;
    }

protected: 

    static bool isValidVector(const vector& v) { return magSqr(v) > VSMALL;}

    
};



}



#endif /* PROPELLER_PROPELLER_SOURCE_ROTOR_GEOMETRY */
