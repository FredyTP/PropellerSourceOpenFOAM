#include "bladeGrid.H"
#include "unitConversion.H"
#include "axisAngleRotation.H"
#include "bladeCell.H"
#include "geometry.H"
namespace Foam
{

defineTypeNameAndDebug(bladeGrid,0);
addToRunTimeSelectionTable(rotorGrid, bladeGrid, dictionary);

bladeGrid::bladeGrid(const dictionary &dict, const rotorGeometry &geometry, const rotorFvMeshSel &rotorFvMeshSel, const bladeModelS* bladeModel, scalar nBlades)
:
    rotorGrid(dict,geometry,rotorFvMeshSel,nBlades)
{
 
    nRadius_ = dict.get<label>("nRadial");
    ijkAddressing::reset(nBlades_,nRadius_,1);
    nChord_ = 1;
    theta_.resize(nBlades_);

    cells_.resize(this->size());
    scalar chord=0;
    bool hasChord = dict.readIfPresent("chord",chord);
    if(hasChord)
    {
        buildBladesConstantChord(chord);
    }
    else if(bladeModel)
    {
        buildBladesFromBladeModel(*bladeModel);
    }
    else
    {
        //Already know that not present, but to generate error msg
        chord = dict.get<scalar>("chord");
    }

    setRotation(getInitialAzimuth());

}

void bladeGrid::assignFvCells()
{
    const auto& cellCenter = meshSel_.mesh().C();
    const auto& weights = meshSel_.mesh().V();
    const auto& cellis = meshSel_.cells();

    const auto& carCS = rotorGeometry_.cartesianCS();
    forAll(cellis,i)
    {
        label celli = cellis[i];
        vector center = carCS.localPosition(cellCenter[celli]);
        List<label> fakeCell({0,1,2,3});
        forAll(cells_,j)
        {
            bladeCell* bCell = dynamic_cast<bladeCell*>(cells_.get(j));
            if(bCell == nullptr)
            {
                continue;
            }       
            if(util::geometry::isInsideCell(bCell->actualPoints(),fakeCell,center))
            {
                cells_.get(j)->addCelli(celli,weights[celli]);          
                break;
            }
        }
    }
}
void bladeGrid::writePythonPlotter(word outputName)
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
        bladeCell* cell = dynamic_cast<bladeCell*>(cells_.get(i));

        const auto& points = cell->actualPoints();
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
void bladeGrid::setRotation(const List<scalar>& thetas)
{
    updateThetas(thetas);
    rotateBlades();
    assignFvCells();
    build();

}
void bladeGrid::buildBladesConstantChord(scalar chord)
{
    List<scalar> radius(nRadius_+1);
    List<scalar> chords(nChord_+1);

    scalar dr = (maxRadius_-minRadius_)/nRadius_;
    scalar dc = (chord)/nChord_;

    //Create radius points
    for(label i = 0; i <nRadius_+1;i++)
    {
        radius[i]=minRadius_+dr*i;
    }

    for(label i = 0; i <nChord_+1;i++)
    {
        chords[i]=-chord/2 + dc*i;
    }

    cells_.resize(ijkAddressing::size());

    for(label ib = 0; ib < nBlades_ ; ib++ )
    {
        for(label ir=0; ir < nRadius_; ir++)
        {
            for(label ic = 0; ic < nChord_; ic++)
            {
                bladeCell* newcell = new bladeCell(rotorGeometry_,radius[ir],radius[ir+1],chords[ic],chords[ic+1]);
                cells_.set(index(ib,ir,ic),newcell);
            }
        }
    }
}
void bladeGrid::buildBladesFromBladeModel(const bladeModelS &bladeModel)
{
    List<scalar> radius(nRadius_+1);

    scalar dr = (maxRadius_-minRadius_)/nRadius_;


    //Create radius points
    for(label i = 0; i <nRadius_+1;i++)
    {
        radius[i]=minRadius_+dr*i;
    }


    cells_.resize(this->size());
    for(label ib = 0; ib < nBlades_ ;ib ++ )
    {
        for(label ir=0; ir<nRadius_; ir++)
        {
            scalar radius0 = radius[ir];
            scalar radius1 = radius[ir+1];
            scalar chord0=0;
            scalar chord1=0;
            scalar twist0=0;
            scalar twist1=0;
            scalar sweep0=0;
            scalar sweep1=0;

            bladeModel.geometryAtRadius(radius0/maxRadius_,chord0,twist0,sweep0);
            bladeModel.geometryAtRadius(radius1/maxRadius_,chord1,twist1,sweep1);

            List<vector> points(4);
            points[0] = vector(radius0,-chord0/2,0);
            points[1] = vector(radius1,-chord1/2,0);
            points[2] = vector(radius1,chord1/2,0);
            points[3] = vector(radius0,chord0/2,0);

            bladeCell* newcell = new bladeCell(rotorGeometry_,points);
            cells_.set(index(ib,ir,0),newcell);

        }
    }

}
void bladeGrid::checkDistribution()
{
    List<scalar> values = chordWiseDistribution(nChord_);
    scalar sum=0;
    scalar tol = 1e-6;
    forAll(values,i)
    {
        sum+=values[i];
    }

    if(std::abs(sum - 1.0) > tol)
    {
        FatalErrorInFunction<<"Chord wise distribution doesnt add up to 1.0"
        <<exit(FatalError);
    }
}
std::function<List<scalar>(label)> bladeGrid::constantDistribution()
{
    return [](label nChord)
    {
        return List<scalar>(nChord,1.0/nChord);
    };
}
void bladeGrid::updateTheta(scalar theta0)
{
    scalar dtheta = constant::mathematical::twoPi/nBlades_;
    forAll(theta_,i)
    {
        theta_[i]=theta0 + dtheta*i;

        if(theta_[i]>constant::mathematical::pi)
        {
            theta_[i]-=constant::mathematical::twoPi;
        }
    }
}

void bladeGrid::updateThetas(const List<scalar> &thetas)
{
    forAll(theta_,i)
    {
        theta_[i]=thetas[i];
        if(theta_[i]>constant::mathematical::pi)
        {
            theta_[i]-=constant::mathematical::twoPi;
        }
    }
}

void bladeGrid::rotateBlades()
{
    for(label ib = 0; ib < nBlades_ ;ib ++ )
    {
        tensor rotation = coordinateRotations::axisAngle::rotation(vector::components::Z,theta_[ib],false);
        for(label ir=0; ir<nRadius_; ir++)
        {
            for(label ic = 0; ic < nChord_; ic++)
            {
                bladeCell* bCell = dynamic_cast<bladeCell*>(cells_.get(index(ib,ir,ic)));
                if(bCell != nullptr)
                {
                    bCell->setRotation(rotation);
                }
            }
        }
    }
}
label bladeGrid::ijkIndex(label iBlade, label iRadius, label iChord)
{
    return ijkAddressing::index(iBlade,iRadius,iChord);
}

List<scalar> bladeGrid::getInitialAzimuth() const
{
    List<scalar> initialAzimuth(nBlades_);
    scalar dtheta = constant::mathematical::twoPi/nBlades_;
    forAll(initialAzimuth,i)
    {
        initialAzimuth[i]=dtheta*i;

        if(initialAzimuth[i]>constant::mathematical::pi)
        {
            initialAzimuth[i]-=constant::mathematical::twoPi;
        }
    }
    return initialAzimuth;
}
}
