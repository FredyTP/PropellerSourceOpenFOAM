#ifndef PROPELLER_ROTOR_MESH_ROTOR_GEOMETRY
#define PROPELLER_ROTOR_MESH_ROTOR_GEOMETRY
#include "scalar.H"
#include "vector.H"

#define NO_RADIUS -1.0

namespace Foam
{
class rotorGeometry
{
private:
    scalar radius_ {NO_RADIUS};
    vector center_ {0,0,0};
    vector direction_ {0,0,0};
    vector psiRef_ {0,0,1}; // not used yet

    bool radiusSet = false;
    bool centerSet = false;
    bool directionSet = false;
    bool psiRefSet = false;

public:

    scalar radius() const {return radius_;}
    const vector& center() const {return center_;}
    const vector& direction() const {return direction_;}
    const vector& psiRef() const {return psiRef_;}

    inline bool isRadiusSet() const { return radiusSet;}
    inline bool isCenterSet() const { return centerSet;}
    inline bool isDirectionSet() const {return directionSet;}
    inline bool isPsiRefSet() const { return psiRefSet;}
   
    /**
     * Read geometry presented data from dictionary
     * If a entry is preset, set the data and set the set flag to true
     * If a entry not present, dont modify the flag
    */
    bool readIfPresent(const dictionary& dict)
    {
        bool ok;

        ok = dict.readIfPresent("radius", radius_);
        radiusSet |= ok;

        ok = dict.readIfPresent("direction",direction_);
        directionSet |= ok;

        ok = dict.readIfPresent("center",center_);
        centerSet |= ok;

        ok = dict.readIfPresent("psiRef",psiRef_);
        psiRefSet |= ok;

        return isReady();
    }
    void setRadius(scalar r)
    {
        radiusSet = true;
        radius_ = r;
    }

    void setCenter(vector c)
    {
        centerSet = true;
        center_ = c;
    }
    void setDirection(vector d)
    {
        directionSet = true;
        direction_ = d;
    }
    void setPsiRef(vector r)
    {
        psiRefSet=true;
        psiRef_ = r;
    }

    inline bool isReady()
    {
        return isSet() && isValid();
    }

    inline bool isValid()
    {
        return (radius_ != NO_RADIUS) &&
               (isValidVector(direction_)) &&
               (isValidVector(psiRef_));
    }

    inline bool isSet()
    {
        return radiusSet && centerSet && directionSet && psiRefSet;
    }

    friend Ostream& operator<<(Ostream& os, const rotorGeometry& r)
    {
        os << "Rotor geometry: "<<endl;
        os << tab <<"- radius: "<<r.radius_<<endl;
        os << tab <<"- center: "<<r.center_<<endl;
        os << tab <<"- direction: "<<r.direction_<<endl;
        os << tab <<"- psiRef: "<< r.psiRef_<<endl;
        return os;
    }
protected: 

    bool isValidVector(vector v){ return magSqr(v) > VSMALL;}

    
};



}



#endif /* PROPELLER_ROTOR_MESH_ROTOR_GEOMETRY */
