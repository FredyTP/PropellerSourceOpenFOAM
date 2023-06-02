#include "polarCell.H"


namespace Foam
{
polarCell::polarCell(const rotorGeometry& rotorGeometry, scalar radius0,scalar radius1, scalar theta0, scalar theta1,label nBlades)
    : gridCell(rotorGeometry), radius0_(radius0), radius1_(radius1), theta0_(theta0), theta1_(theta1)
{
    dr_ = this->radius1()-this->radius0();
    dt_ = this->theta1()-this->theta0();

    setCenter(getCellCenter());

    factor_ = nBlades * this->dr() * this->dt() / constant::mathematical::twoPi;
}
vector polarCell::getCellCenter() const
{
    return vector(this->radius0()+dr_/2,this->theta0()+dt_/2,0);
}

}
