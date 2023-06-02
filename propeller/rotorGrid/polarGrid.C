#include "polarGrid.H"
#include "mathematicalConstants.H"
#include "regularInterpolation.H"
#include "polarCell.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{


defineTypeNameAndDebug(polarGrid,0);
addToRunTimeSelectionTable(rotorGrid, polarGrid, dictionary);

polarGrid::polarGrid(const dictionary &dict, const rotorGeometry &geometry, const rotorFvMeshSel &rotorFvMeshSel, const bladeModelS* bladeModel, scalar nBlades)
: rotorGrid(dict, geometry, rotorFvMeshSel, nBlades)
{
    label nRadius = dict.get<label>("nRadial");
    label nTheta = dict.get<label>("nAzimutal");
    radius_.resize(nRadius+1);
    theta_.resize(nTheta+1);
    ijkAddressing::reset(nRadius,nTheta,0);

    scalar dr = (maxRadius_-minRadius_)/nRadius;
    scalar dt = (constant::mathematical::twoPi)/nTheta;


    for(label i = 0; i <nRadius+1;i++)
    {
        radius_[i]=minRadius_+dr*i;
    }


    for(label i = 0; i <nTheta+1;i++)
    {
        theta_[i]=-constant::mathematical::pi+dt*i;
    }

    buildGrid();
    assignFvCells();
    build();

}

void polarGrid::assignFvCells()
{
    const auto& cellCenter = meshSel_.mesh().C();
    const auto& weights = meshSel_.mesh().V();
    const auto& cellis = meshSel_.cells();

    const auto& cylCS = rotorGeometry_.cylindricalCS();
    forAll(cellis,i)
    {
        label celli = cellis[i];
        vector cc = cellCenter[celli];
        vector polar = cylCS.localPosition(cc);
        label ir=0;
        label unused;
        label it=0;
        int result = regularInterpolation<scalar,scalar,1>::FindIndex(polar.x(),this->radius(),ir,unused);
        regularInterpolation<scalar,scalar,1>::FindIndex(polar.y(),this->theta(),it,unused);

        if(result == 1)
        {
            this->cell(ir,it).addCelli(celli,weights[celli]);
        }
    }
}

void polarGrid::build()
{
    forAll(cells_,i)
    {
        cells_[i].checkCells();
        cells_[i].buildWeigths();
    }
    updateCenters();
}


void polarGrid::buildGrid()
{
    label radialCells = radius_.size()-1;
    label thetaCells = theta_.size()-1;

    cells_.resize(radialCells*thetaCells);

    for(label i = 0; i< radialCells;i++)
    {
        for(label j = 0; j<thetaCells;j++)
        {
            polarCell* newcell = new polarCell(rotorGeometry_,radius_[i],radius_[i+1],theta_[j],theta_[j+1],nBlades_);
            cells_.set(index(i,j,0),newcell);
        }
    }
}

}

