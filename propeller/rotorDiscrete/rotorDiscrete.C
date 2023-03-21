#include "rotorDiscrete.H"
#include "cellSet.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "unitConversion.H"
#include "rotorGeometry.H"
#include "delaunayTriangulation.H"
#include <fstream>


namespace Foam
{

defineTypeNameAndDebug(rotorDiscrete,0);



rotorDiscrete::rotorDiscrete()
{

}
void rotorDiscrete::buildCoordinateSystem(const rotorGeometry& geometry)
{
    rotorGeometry_ = geometry;

    cylCS_ = coordSystem::cylindrical
                (
                    rotorGeometry_.center, //centerd to local
                    rotorGeometry_.direction, //z-axis 
                    rotorGeometry_.psiRef  //x-axis
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
void rotorDiscrete::fromRotorMesh(const rotorMesh &rotorMesh)
{
    Info<< "Building rotor Discrete from mesh" <<endl;
    discreteMode_ = discreteMode::dmMesh;
    const labelList cells = rotorMesh.cells();
    const auto &mesh = rotorMesh.mesh();
    
    cylPoints_.resize(cells.size());
    localBlade_.resize(cells.size());


    volScalarField volume
    (
        IOobject
        (
            "propeller:rotorVolume",
            rotorMesh.mesh().time().timeName(),
            rotorMesh.mesh()
        ),
        rotorMesh.mesh(),
        dimensionedScalar(dimVolume, Zero)
    );

    forAll(cells, i)
    {
        label celli = cells[i];
        vector cellCenter = mesh.C()[celli];

        //Project points over the disk
        //From global to local cyl position
        // Global -> local cyl
        cylPoints_[i] = cylCS_.localPosition(cellCenter);
        cylPoints_[i].z()=0;
        
        //Not the most efficient way to compute global coords
        localBlade_[i] = this->bladeLocalFromPoint(cylPoints_[i]);
        volume[celli] = mesh.V()[celli];
    }

    volume.write();



    //-------------BEGIN JUST TEST-------------//

        /*List<XYZ> ps(6);
        ps[0] = XYZ{-0.1,-0.1,0};
        ps[1] = XYZ{0.1,-0.1,0};
        ps[2] = XYZ{-0.1,0.1,0};
        ps[3] = XYZ{0.1,0.1,0};
        ps[4] = XYZ{0,0,0};
        ps[5] = XYZ{0.15,0,0};*/
        List<point> ps(cylPoints_.size());
        for(label i = 0; i <cylPoints_.size(); i++)
        {
            point pt = cylCS_.globalPosition(cylPoints_[i]);
            ps[i] = point{pt.x(),pt.z(),0};
        }

        std::sort(ps.begin(),ps.end(),[](const point& p1, const point& p2)
            {
                return p1.x() < p2.x();
            });

        List<point> orto;
        List<List<label>> voronoid;
        delaunayTriangulation::Voronoid
        (
            ps,
            orto,
            voronoid,
            delaunayTriangulation::circularRegion(0.15),
            delaunayTriangulation::intersectCircle(0.15)
        );
        
        std::string x_string = "x = [";
        std::string y_string = "y = [";
        std::string xc_string = "xc = [";
        std::string yc_string = "yc = [";

        for(label i = 0 ; i < ps.size();i++)
        {
            xc_string += std::to_string(ps[i].x());
            
            yc_string += std::to_string(ps[i].y());

            if(i != ps.size()-1)
            {
                xc_string +=",";
                yc_string +=",";
            }
        }

        xc_string +="]";
        yc_string +="]";

        for(label i = 0 ; i < orto.size();i++)
        {
            x_string += std::to_string(orto[i].x());
            
            y_string += std::to_string(orto[i].y());

            if(i != orto.size()-1)
            {
                x_string +=",";
                y_string +=",";
            }
        }

        x_string +="]";
        y_string +="]";

        std::string tri_str = "tri = [";
        for(label i = 0; i< voronoid.size(); i++)
        {
            auto vor = voronoid[i];
            if(vor.size()==0) continue;
            tri_str += "[";   
            for(label j=0; j <vor.size();j++)
            {
                tri_str += std::to_string(vor[j]);
                tri_str += ",";
            }
            tri_str += std::to_string(vor[0]);
            tri_str += "]";
            if(i!= voronoid.size()-1)
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
    /*wchar_t *program = Py_DecodeLocale("./", NULL);
    if (program == NULL) {
        fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
        //exit(1);
    }
    Py_SetProgramName(program); */ /* optional but recommended */
    /*Py_Initialize();
    PyRun_SimpleString(pyplot.c_str());
    if (Py_FinalizeEx() < 0) {
        //exit(120);
    }
    PyMem_RawFree(program);*/
    


     //-------------END JUST TEST-------------//
}
bool rotorDiscrete::read(const dictionary &dict)
{
    return false;
}

}