#include "gridCell.H"


namespace Foam
{
gridCell::gridCell(label iRadius, label iTheta, const List<scalar> &radius, const List<scalar> &theta)
    : iRadius_(iRadius), iTheta_(iTheta), radius_(radius), theta_(theta)
{
    dr_ = radius1()-radius0();
    dt_ = theta1()-theta0();

    center_ = vector(radius0()+dr_/2,theta0()+dt_/2,0);
}

void gridCell::addCelli(label celli,scalar area)
{
    cellis_.append(celli);
    weights_.append(area);
}
void gridCell::build()
{
    scalar totalArea=0.0;
    if(cellis_.size()==0)
    {
        //FatalErrorInFunction
         //   <<"Some gridCells doesnt contain any fvCell"
         //   <<exit(FatalError);
    }
    
    forAll(weights_,i)
    {
        totalArea+=weights_[i];
    }

    forAll(weights_,i)
    {
        weights_[i]/=totalArea;
    }


}

}
