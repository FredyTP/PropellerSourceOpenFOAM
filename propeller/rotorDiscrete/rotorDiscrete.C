#include "rotorDiscrete.H"
#include "cellSet.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "unitConversion.H"
#include "rotorGeometry.H"
#include "delaunayTriangulation.H"
#include <fstream>
#include "rotorTriCell.H"

namespace Foam
{

defineTypeNameAndDebug(rotorDiscrete,0);



rotorDiscrete::rotorDiscrete()
{

}
void rotorDiscrete::buildCoordinateSystem(const rotorGeometry& geometry)
{
    rotorGeometry_ = geometry;
    carCS_ = coordSystem::cartesian
                (
                    rotorGeometry_.center(), //centerd to local
                    rotorGeometry_.direction(), //z-axis 
                    rotorGeometry_.psiRef()  //x-axis
                );

    cylCS_ = coordSystem::cylindrical //local Cartesian to cylindrical
                (
                    rotorGeometry_.center(), //centerd to local
                    rotorGeometry_.direction(), //z-axis 
                    rotorGeometry_.psiRef()  //x-axis
                );

    carToCylCS_ = coordSystem::cylindrical
                    (
                        vector(0,0,0), //same center
                        vector(0,0,1), //same z-axis
                        vector(1,0,0)  //same x-axis
                    );

}
tensor rotorDiscrete::bladeLocalFromPoint(const point &localPoint) const
{
    //z- up, y -outwards from center, x perpendicular y,z (leading edge to trailing edge)
    point global,origin;

    global = cylCS_.globalPosition(localPoint);
    origin = cylCS_.origin();
    
    tensor rotTensor(cylCS_.R());

    //z-axis
    const vector ax3 = rotTensor.col<2>(); // == e3 (already normalized)

    //y-axis (radial direction)
    vector ax2(global - origin);

    ax2.removeCollinear(ax3);

    const scalar magAxis2(mag(ax2));

    // Trap zero size and colinearity
    if (magAxis2 < SMALL)
    {
        return rotTensor;
    }

    ax2 /= magAxis2;  // normalise

    // Replace with updated local axes

    rotTensor.col<0>(ax2^ax3);
    rotTensor.col<1>(ax2);

    return rotTensor;
}
void rotorDiscrete::fromRotorMesh(const rotorFvMeshSel &rotorFvMeshSel, word integration)
{
    Info<< "Building rotor Discrete from mesh" <<endl;
    discreteMode_ = discreteMode::dmMesh;

    //Selected rotor radius (real used, no from mesh)
    const scalar radius = rotorGeometry_.radius();
    const scalar sqrRadius = radius * radius;
    const scalar idealArea = constant::mathematical::pi * sqrRadius;

    //List ref.
    const labelList& cells = rotorFvMeshSel.cells(); 
    const vectorField& cellCenter = rotorFvMeshSel.mesh().C();

    List<point> centers(cells.size());
    List<label> cellis(cells.size());
    label iCent = 0;
    forAll(cells,i)
    {
        label celli = cells[i]; //cell index

        //Get local center coordinates
        point testCenter = carCS_.localPosition(cellCenter[celli]);
        testCenter.z() = 0; //proyect over plane
        if(magSqr(testCenter) <= sqrRadius)
        {
            centers[iCent] = testCenter;
            cellis[iCent] = celli;
            iCent++;
        }

    }

    //Resize to fit correct data
    centers.resize(iCent);
    cellis.resize(iCent);

    List<List<label>> voroCells;
    List<point> vertex;

    //Create voronoid diagram of proyected centers for 2D meshing
    delaunayTriangulation::Voronoid
    (
        centers,
        vertex,  
        voroCells, 
        delaunayTriangulation::circularRegion(radius),
        delaunayTriangulation::intersectCircle(radius)
    );

    //Add vertex points
    carPoints_ = vertex;
    carPoints_.resize(vertex.size() + centers.size());
    integrationPoints_.resize(carPoints_.size());
    rotorCells_.resize(centers.size());

    area_ = 0.0;
    scalar ta = 0.0;
    // add cell centers to rotor points and create cells
    forAll(voroCells,i)
    {
        label celli = cellis[i]; //Now indexing is from cellis (used cells)
        carPoints_[i+vertex.size()] = centers[i];
        rotorCells_.set(i,rotorCell::New(integration,i+vertex.size(),voroCells[i],carPoints_,celli));
        rotorCells_[i].updateIntegrationList(integrationPoints_);
        area_ += rotorCells_[i].area();
        ta += (static_cast<rotorTriCell*>(rotorCells_.get(i)))->areaTri();
    }   

    Info<< "Total disk area: "<<area_<<endl;
    Info<< "Total diak tri area: "<<ta<<endl;
    Info<< "Idea disk area: "<<idealArea<<endl;

    //TODO: ?Add aditional integration points
    
    //Add cyl coord system and local blade tensor

    cylPoints_.resize(carPoints_.size());
    localBlade_.resize(carPoints_.size());
    forAll(carPoints_,i)
    {
        cylPoints_[i] = carToCylCS_.localPosition(carPoints_[i]);
        localBlade_[i] = this->bladeLocalFromPoint(cylPoints_[i]);
    }


//-------------BEGIN JUST TEST-------------//


        std::string x_string = "x = [";
        std::string y_string = "y = [";
        std::string xc_string = "xc = [";
        std::string yc_string = "yc = [";

        for(label i = 0 ; i < centers.size();i++)
        {
            xc_string += std::to_string(centers[i].x());
            
            yc_string += std::to_string(centers[i].y());

            if(i != centers.size()-1)
            {
                xc_string +=",";
                yc_string +=",";
            }
        }

        xc_string +="]";
        yc_string +="]";

        for(label i = 0 ; i < vertex.size();i++)
        {
            x_string += std::to_string(vertex[i].x());
            
            y_string += std::to_string(vertex[i].y());

            if(i != vertex.size()-1)
            {
                x_string +=",";
                y_string +=",";
            }
        }

        x_string +="]";
        y_string +="]";

        std::string tri_str = "tri = [";
        for(label i = 0; i< voroCells.size(); i++)
        {
            auto vor = voroCells[i];
            if(vor.size()==0) continue;
            tri_str += "[";   
            for(label j=0; j <vor.size();j++)
            {
                tri_str += std::to_string(vor[j]);
                tri_str += ",";
            }
            tri_str += std::to_string(vor[0]);
            tri_str += "]";
            if(i!= voroCells.size()-1)
            {
                tri_str += ",";
            }

        }
        tri_str += "]";

    std::string pyplot= "import numpy as np;\n";
    pyplot+= "import matplotlib.pyplot as plot;\n";
    pyplot += "import matplotlib as mp;\n";
    pyplot+= xc_string + "\n";
    pyplot+= yc_string + "\n";
    pyplot+= x_string + "\n";
    pyplot+= y_string +"\n";
    pyplot+= tri_str + "\n";

    pyplot+="tris = []\n";
    pyplot+="for i in range(len(tri)):\n";
    pyplot+="\txp = np.zeros(len(tri[i]))\n";
    pyplot+="\typ = np.zeros(len(tri[i]))\n";
    pyplot+="\tfor j in range(len(tri[i])):\n";
    pyplot+="\t\txp[j]=x[tri[i][j]]\n";
    pyplot+="\t\typ[j]=y[tri[i][j]]\n";
    //pyplot+="\tif(len(tri[i])>2):\n";
    pyplot+="\tplot.plot(xp,yp)\n";
    pyplot+="\n";
    pyplot+="plot.plot(xc,yc,marker='o',linewidth = 0)\n";
    pyplot+="plot.show()\n";
    std::ofstream file("triangulation.py",std::ios::out);
    file<<pyplot;
    file.close();

     //-------------END JUST TEST-------------//
}
bool rotorDiscrete::read(const dictionary &dict)
{
    return false;
}

}