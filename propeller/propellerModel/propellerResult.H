#ifndef PROPELLER_PROPELLER_MODEL_PROPELLER_RESULT
#define PROPELLER_PROPELLER_MODEL_PROPELLER_RESULT

#include "fvMesh.H"
#include "regIOobject.H"
namespace Foam
{

/**
 * @brief result integrated data information
 * Adimensional parameters defined as website : 
 * https://m-selig.ae.illinois.edu/props/propDB.html
 * suggests
*/
class propellerResult
{   
public:
    vector force {0,0,0}; //- x: y: z: thrust in local coordinates
    vector torque {0,0,0};    //- Aerodynamic torque on rotor axis
    scalar power {0};     //- Power needed to rotate the propeller
    scalar J{0};          //- Adimensional: refV / (n*D)
    scalar eta{0};        //- Thrust* refV / power
    scalar CT{0};         //- Thrust coefficient T/(rho*n^2*D^4)
    scalar CQ{0};         //- Torque coefficient Q/(rho*n^2*D^5)
    scalar CP{0};         //- Power coefficient Power/(rho*n^3*D^5)

    propellerResult()
    {
        
    }
    /**
     * @brief Updates internal data
     * 
     * @param omega 
     * @param speedRef 
     * @param rhoRef 
     * @param radius 
     */
    void update(scalar omega, scalar speedRef,scalar rhoRef, scalar radius)
    {
        power = torque.z() * omega;
        updateEta(speedRef);
        updateJ(speedRef,omega,radius);
        updateCT(rhoRef,omega,radius);
        updateCP(rhoRef,omega,radius);
    }
    inline scalar thrust() const {return force.z();}

    static Ostream& OutputDefinition(Ostream& os)
    {
        if(!definitionShown_)
        {
            definitionShown_ = true;
            indent(os) << "Rotor adimiensional parameters definitions: "<<endl;
            os.incrIndent();
            indent(os) <<"- J (Advance Coefficient):  Vref/(n*D)"<<endl;
            indent(os) <<"- Eta (Efficiency): (T * Vref)/P"<<endl;
            indent(os) <<"- CT (Thrust Coefficient): T/(rho * n^2 * D^4)"<<endl;
            indent(os) <<"- CQ (Torque Coefficient): Q/(rho * n^2 * D^5)"<<endl;
            indent(os) <<"- CP (Power Coefficient): P/(rho * n^3 * D^5)"<<endl;
            os.decrIndent();
            indent(os) << "Where: "<<endl;
            os.incrIndent();
            indent(os)<< "- Vref is reference velocity"<<endl;
            indent(os)<< "- n is rotational speed in revolution/second"<<endl;
            indent(os)<< "- D is rotor diameter"<<endl;
            indent(os)<< "- T is rotor thrust"<<endl;
            indent(os)<< "- Q is rotor torque"<<endl;
            indent(os)<< "- P is rotor power"<<endl;
            indent(os)<< "- rho is reference density"<<endl;
            os.decrIndent();
        }
        return os;
    }

    friend Ostream& operator<<(Ostream& os, const propellerResult& pr)
    {   
        os.incrIndent();
        indent(os) <<"- Force = "<<pr.force<<endl;
        indent(os) <<"- Torque = "<<pr.torque<<endl;
        indent(os) <<"- Power = "<<pr.power<<endl;
        indent(os) <<"- J = "<<pr.J<<endl;
        indent(os) <<"- Eta = "<<pr.eta<<endl;
        indent(os) <<"- CT = "<<pr.CT<<endl;
        indent(os) <<"- CQ = "<<pr.CQ<<endl;
        indent(os) <<"- CP = "<<pr.CP<<endl; 
        os.decrIndent();
        return os;
    }

    /**
     * Update Js with force and power already set
     * J defined as Adimensional: refV / (n*D)
    */
    void updateJ(scalar vRef, scalar omega, scalar radius)
    {
        scalar vRot = omega/constant::mathematical::pi * radius;
        if(vRot == 0)
        {
            J = 0;
        }
        else
        {
            J = vRef/vRot;
        }
    }

    /**
     * @brief: Update eta value with force and power already set
     * eta defined as Adimensional: Thrust * vRef / power
    */
    void updateEta(scalar vRef)
    {
        if(power == 0)
        {
            eta = 0;
        }
        else
        {
            eta=thrust() * vRef / power;
        }
    }

    void updateCT(scalar rhoRef,scalar omega, scalar radius)
    {
        scalar den = rhoRef * 
                pow((omega / constant::mathematical::twoPi),2)*
                pow(2*radius,4);
        if(den == 0)
        {
            CT=0;
        }
        else
        {
            CT = thrust()/den;
        }
    }

    void updateCP(scalar rhoRef, scalar omega,scalar radius)
    {
        scalar den = rhoRef * 
                pow((omega / constant::mathematical::twoPi),3)*
                pow(2*radius,5);
        if(den == 0)
        {
            CP=0;
        }
        else
        {
            CP = power/den;
        }
        CQ = CP/constant::mathematical::twoPi;
    }

private:
    static bool definitionShown_;
};

/**
 * @brief Class enabling auto input-output of propellerResult
 * 
 */
class regIOpropellerResult
    : public propellerResult, public regIOobject
{

public:
    regIOpropellerResult(word name,const fvMesh& mesh)
        : regIOobject
        (
            IOobject
            (
                name+":result",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            )
        )
    {

    }
    regIOpropellerResult& operator=(propellerResult result)
    {
        this->force = result.force;
        this->torque = result.torque;
        this->power = result.power;
        this->J = result.J;
        this->eta = result.eta;
        this->CT = result.CT;
        this->CQ = result.CQ;
        this->CP = result.CP;

        return *this;
    }
    bool writeData(Ostream& os) const
    {
        
        os <<"force = "<<this->force<<";"<<endl;
        os <<"thrust = "<<this->thrust()<<";"<<endl;
        os <<"torque = "<<this->torque<<";"<<endl;
        os <<"power = "<<this->power<<";"<<endl;
        os <<"j = "<<this->J<<";"<<endl;
        os <<"eta = "<<this->eta<<";"<<endl;
        os <<"ct = "<<this->CT<<";"<<endl;
        os <<"cq = "<<this->CQ<<";"<<endl;
        os <<"cp = "<<this->CP<<";"<<endl; 
        
        return true;
    }



};

}

#endif /* PROPELLER_PROPELLER_MODEL_PROPELLER_RESULT */
