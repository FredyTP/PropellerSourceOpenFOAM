#include "rotorDiscrete.H"
#include "cellSet.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "unitConversion.H"
#include "rotorGeometry.H"
#include "delaunayTriangulation.H"
#include <fstream>
#include "syncTools.H"
#include "rotorGrid.H"
#include "regularInterpolation.H"

namespace Foam
{

defineTypeNameAndDebug(rotorDiscrete, 0);

const Enum
<
    rotorDiscrete::discreteMethod
>
rotorDiscrete::discreteMethodNames_
({
    {discreteMethod::dmeVoronoid, "voronoid"},
    {discreteMethod::dmeIntersection, "intersection"}
});



rotorDiscrete::rotorDiscrete(const dictionary& dict)
{
    this->read(dict);
}
void rotorDiscrete::createGrid(const rotorGeometry &geometry)
{
    rotorGeometry_ = &geometry;
    
    carCS_ = coordSystem::cartesian(
        rotorGeometry_->center().get(),    // centerd to local
        rotorGeometry_->direction().get(), // z-axis
        rotorGeometry_->psiRef().get()     // x-axis
    );

    cylCS_ = coordSystem::cylindrical // local Cartesian to cylindrical
        (
            rotorGeometry_->center().get(),    // centerd to local
            rotorGeometry_->direction().get(), // z-axis
            rotorGeometry_->psiRef().get()     // x-axis
        );

    carToCylCS_ = coordSystem::cylindrical(
        vector(0, 0, 0), // same center
        vector(0, 0, 1), // same z-axis
        vector(1, 0, 0)  // same x-axis
    );

    grid_ = rotorGrid
    (
        nRadial,
        nAzimutal,
        rotorGeometry_->innerRadius().get(),
        rotorGeometry_->radius().get(),
        cylCS_
    );
}


void rotorDiscrete::writePythonPlotter(word process)
{
    /*std::string x_string = "x = [";
    std::string y_string = "y = [";
    std::string xc_string = "xc = [";
    std::string yc_string = "yc = [";

    for (label i = 0; i < grid.cells().size() i++)
    {
        const auto& gridcell = grid.cells()[i];

        vector center = gridcell.center();
        center = carToCylCS_.globalPosition(center);
        xc_string += std::to_string(center.x());

        yc_string += std::to_string(center.y());

        if (i != rotorCells_.size() - 1)
        {
            xc_string += ",";
            yc_string += ",";
        }
    }

    xc_string += "]";
    yc_string += "]";


    std::string trix_str = "trix = [";
    std::string triy_str = "triy = [";
    for (label i = 0; i < grid.cells().size(); i++)
    {
        const auto& gridcell = grid.cells()[i];
        if (vor.size() == 0)
            continue;
        tri_str += "[";
            tri_str += std::to_string(vor[j]);
            tri_str += ",";

        tri_str += std::to_string(vor[0]);
        tri_str += "]";
        if (i != rotorCells_.size() - 1)
        {
            tri_str += ",";
        }
    }
    tri_str += "]";

    std::string pyplot = "import numpy as np;\n";
    pyplot += "import matplotlib.pyplot as plot;\n";
    pyplot += "import matplotlib as mp;\n";
    pyplot += xc_string + "\n";
    pyplot += yc_string + "\n";
    pyplot += x_string + "\n";
    pyplot += y_string + "\n";
    pyplot += tri_str + "\n";

    pyplot += "tris = []\n";
    pyplot += "for i in range(len(trix)):\n";
    pyplot += "\txp = np.zeros(len(trix[i]))\n";
    pyplot += "\typ = np.zeros(len(trix[i]))\n";
    pyplot += "\tfor j in range(len(trix[i])):\n";
    pyplot += "\t\txp[j]=trix[i][j]\n";
    pyplot += "\t\typ[j]=triy[i][j]\n";
    // pyplot+="\tif(len(tri[i])>2):\n";
    pyplot += "\tplot.plot(xp,yp)\n";
    pyplot += "\n";
    pyplot += "plot.plot(xc,yc,marker='o',linewidth = 0)\n";
    pyplot += "plot.show()\n";
    std::ofstream file("rotorGrid" + std::to_string(Pstream::myProcNo()) + ".py", std::ios::out);
    file << pyplot;
    file.close();*/
}
void rotorDiscrete::setFvMeshSel(const rotorFvMeshSel &rotorFvMeshSel)
{
    Info<<endl;
    Info << "Assigning fvMeshCell to rotorGrid: " << endl;
    rotorMeshSel_ = &rotorFvMeshSel;
    // List ref.
    const vectorField& cellCenter = rotorFvMeshSel.mesh().C();
    const scalarField& cellVol = rotorFvMeshSel.mesh().V();
    const labelList& cellis = rotorFvMeshSel.cells();

    grid_.assignFvCells(cellCenter,cellVol,cellis);
    grid_.setCenterFromClosestCell(cellCenter);
    grid_.build();
}

bool rotorDiscrete::read(const dictionary &dict)
{

    //integrationMode_ = rotorCell::integrationModeNames_.getOrDefault("integrationMode",dict,rotorCell::integrationMode::imCenter);
    //discreteMethod_ = discreteMethodNames_.getOrDefault("discreteMethod",dict,discreteMethod::dmeVoronoid);

    //For voronoid method vertex cells are not needed, but for intersection is recommended
    //bool defaultVertex = discreteMethod_ == discreteMethod::dmeIntersection;
    //includeVertex_ = dict.getOrDefault<bool>("includeIfVertex",defaultVertex);

    //refinementLevel_ = dict.getOrDefault<label>("borderRefinement",0);

    //correctCenters_ = dict.getOrDefault<bool>("correctCenters",false);
    
    nRadial = dict.get<label>("nRadial");
    nAzimutal = dict.get<label>("nAzimutal");

    Info<<endl;    
    Info << "Reading rotor Discrete dict:" << endl;
    Info.stream().incrIndent();
    indent(Info)<<"- Radial Cells: "<<nRadial<<endl;
    indent(Info)<<"- Azimutal Cells: "<<nAzimutal<<endl;

    //indent(Info)<<"- Discrete method: "<<rotorDiscrete::discreteMethodNames_.get(discreteMethod_)<<endl;
    //indent(Info)<<"- Border refinement: "<<refinementLevel_<<endl;
    //indent(Info)<<"- Correct centers: "<<correctCenters_<<endl;
    //indent(Info)<<"- Integration mode: "<<rotorCell::integrationModeNames_.get(integrationMode_)<<endl;
    //indent(Info)<<"- Include if vertex: "<<includeVertex_<<endl;
    Info.stream().decrIndent();

    return true;
}

}
