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
#include "bladeGrid.H"
#include "regularInterpolation.H"

namespace Foam
{

defineTypeNameAndDebug(rotorDiscrete, 0);

const Enum
<
    rotorDiscrete::sampleMode
>
rotorDiscrete::sampleModeNames_
({
    {sampleMode::spCenter, "center"},
    {sampleMode::spClosestCell, "closestCell"},
});


rotorDiscrete::rotorDiscrete(const dictionary& dict)
   : dict_(dict)
{
    this->read(dict);
    
}
void rotorDiscrete::setGeometry(const rotorGeometry &geometry, scalar nBlades)
{
    rotorGeometry_ = &geometry;
    nBlades_ = nBlades;
    carCS_ = geometry.cartesianCS();
    cylCS_ = geometry.cylindricalCS();
    carToCylCS_ = geometry.cartesianToCylindrical();
}

void rotorDiscrete::updateTheta(scalar theta)
{
    bladeGrid* bGrid = dynamic_cast<bladeGrid*>(grid_.get());
    if(bGrid)
    {
        bGrid->setRotation(theta);
    }
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
void rotorDiscrete::setFvMeshSel(const rotorFvMeshSel &rotorFvMeshSel, const bladeModelS& bladeModel)
{
    Info<<endl;
    Info << "Assigning fvMeshCell to rotorGrid: " << endl;
    rotorMeshSel_ = &rotorFvMeshSel;

    grid_ = rotorGrid::New(dict_,*rotorGeometry_,rotorFvMeshSel, bladeModel, nBlades_);
    // List ref.

    Info<<"Assigning fv cells"<<endl;

    grid_->assignFvCells();
    if(sampleMode_ == spClosestCell)
    {
        grid_->setCenterFromClosestCell();
    }
    Info<<"Building "<<endl;
    grid_->build();

    volScalarField selected
    (
        IOobject
        (
            "selectedCells",
            rotorMeshSel_->mesh().time().timeName(),
            rotorMeshSel_->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rotorMeshSel_->mesh(),
        dimensionedScalar(dimless, Zero)
    );

    forAll(grid_->cells(),i)
    {
        grid_->cells()[i].applyField<scalar>(selected.ref(false),1);
    }

    selected.write();

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
    sampleMode_ = sampleModeNames_.getOrDefault("sampleMode",dict,sampleMode::spClosestCell);

    Info<<endl;    
    Info << "Reading rotor Discrete dict:" << endl;
    Info.stream().incrIndent();
    indent(Info)<<"- Radial Cells: "<<nRadial<<endl;
    indent(Info)<<"- Azimutal Cells: "<<nAzimutal<<endl;
    indent(Info)<<"- Sample Mode: "<<sampleModeNames_.get(sampleMode_)<<endl;
    //indent(Info)<<"- Discrete method: "<<rotorDiscrete::discreteMethodNames_.get(discreteMethod_)<<endl;
    //indent(Info)<<"- Border refinement: "<<refinementLevel_<<endl;
    //indent(Info)<<"- Correct centers: "<<correctCenters_<<endl;
    //indent(Info)<<"- Integration mode: "<<rotorCell::integrationModeNames_.get(integrationMode_)<<endl;
    //indent(Info)<<"- Include if vertex: "<<includeVertex_<<endl;
    Info.stream().decrIndent();

    return true;
}

}
