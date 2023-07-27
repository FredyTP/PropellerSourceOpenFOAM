#include "polarGrid.H"
#include "mathematicalConstants.H"
#include "polarCell.H"
#include "addToRunTimeSelectionTable.H"
#include "RegularInterpolation.H"

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
    ijkAddressing::reset(nRadius,nTheta,1);

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
void polarGrid::writePythonPlotter(word outputName)
{

    word procName="";
    if(Pstream::parRun())
    {
        procName = std::to_string(Pstream::myProcNo());
    }   
    

    /**VERTEX*/
    std::string x_string = "x = [";
    std::string y_string = "y = [";


    for (label i = 0; i < cells_.size(); i++)
    {
        polarCell* cell = dynamic_cast<polarCell*>(cells_.get(i));

        scalar r0 = cell->radius0();
        scalar r1 = cell->radius1();
        scalar t0 = cell->theta0();
        scalar t1 = cell->theta1();

        vector p1(r0,t0,0);
        vector p2(r0,t1,0);
        vector p3(r1,t1,0);
        vector p4(r1,t0,0);
        List<vector> points;
        points.append(p1);
        label reso=10;
        scalar dt = (t1-t0)/reso;
        for(label i = 0; i < reso;i++)
        {
            p1+=vector(0,dt,0);
            points.append(p1);
        }
        points.append(p2);
        points.append(p3);
        for(label i = 0; i < reso;i++)
        {
            p3-=vector(0,dt,0);
            points.append(p3);
        }
        points.append(p4);

        forAll(points,i)
        {
            points[i]=rotorGeometry_.cartesianToCylindrical().globalPosition(points[i]);
        }
        
        x_string+="[";
        y_string+="[";

        for (label j = 0; j < points.size(); j++)
        {
            x_string += std::to_string(points[j].x());
            y_string += std::to_string(points[j].y());

            x_string += ",";
            y_string += ",";
            
        }
        x_string += std::to_string(points[0].x());
        y_string += std::to_string(points[0].y());

        x_string+="]";
        y_string+="]";

        if (i != cells_.size() - 1)
        {
            x_string += ",";
            y_string += ",";
        }
    }

    x_string += "]";
    y_string += "]";


    /**CENTERS*/
    std::string xc_string = "xc = [";
    std::string yc_string = "yc = [";
    for (label i = 0; i < cells_.size(); i++)
    {

        vector center = cells_[i].center();
        center = rotorGeometry_.cartesianToCylindrical().globalPosition(center);
    
        xc_string += std::to_string(center.x());
        yc_string += std::to_string(center.y());

        if (i != cells_.size() - 1)
        {
            xc_string += ",";
            yc_string += ",";
        }
    }

    xc_string += "]";
    yc_string += "]";

    std::string pyplot = "import numpy as np;\n";
    pyplot += "import matplotlib.pyplot as plot;\n";
    pyplot += "import matplotlib as mp;\n";
    pyplot += xc_string + "\n";
    pyplot += yc_string + "\n";
    pyplot += x_string + "\n";
    pyplot += y_string + "\n";

    pyplot += "for i in range(len(x)):\n";
    pyplot += "\tplot.plot(x[i],y[i],'k')\n";
    pyplot += "\n";
    pyplot += "plot.plot(xc,yc,marker='o',linewidth = 0)\n";
    pyplot += "plot.gca().set_aspect('equal')\n";
    pyplot += "plot.show()\n";
    std::ofstream file(outputName + procName + ".py", std::ios::out);
    file << pyplot;
    file.close();
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
        int result = RegularInterpolation<scalar,scalar,1>::FindIndex(polar.x(),this->radius(),ir,unused);
        RegularInterpolation<scalar,scalar,1>::FindIndex(polar.y(),this->theta(),it,unused);

        if(result == 1)
        {
            this->cell(ir,it).addCelli(celli,weights[celli]);
        }
    }
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

