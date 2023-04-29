#ifndef PROPELLER_ROTOR_DISCRETE_GRID_CELL
#define PROPELLER_ROTOR_DISCRETE_GRID_CELL

#include "scalar.H"
#include "List.H"
#include "vector.H"

namespace Foam
{
class gridCell
{
public:
    gridCell(label iRadius, label iTheta, const List<scalar>& radius, const List<scalar>& theta);


    void addCelli(label celli,scalar area);

    void build();

    const labelList& cellis() const { return cellis_; }
    const List<scalar>& weights() const {return weights_;}

    scalar radius()const {return center_.x();}
    scalar theta() const {return center_.y();}

    const vector& center() const {return center_;}
    scalar dr() const { return dr_;}
    scalar dt() const {return dt_;}

    scalar radius0(){return radius_[iRadius_];}
    scalar radius1(){return radius_[iRadius_+1];}
    scalar theta0(){return theta_[iTheta_];}
    scalar theta1(){return theta_[iTheta_+1];}
    
private:
    vector center_; //r-theta-0

    label iRadius_; //i and i+1
    label iTheta_; //i and i+1

    scalar dr_;
    scalar dt_;

    List<label> cellis_;
    List<scalar> weights_;

    const List<scalar>& radius_;
    const List<scalar>& theta_;
};

}

#endif /* PROPELLER_ROTOR_DISCRETE_GRID_CELL */
