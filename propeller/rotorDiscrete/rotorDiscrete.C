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

defineTypeNameAndDebug(rotorDiscrete, 0);

rotorDiscrete::rotorDiscrete()
{
}
void rotorDiscrete::buildCoordinateSystem(const rotorGeometry &geometry)
{
    rotorGeometry_ = geometry;
    carCS_ = coordSystem::cartesian(
        rotorGeometry_.center(),    // centerd to local
        rotorGeometry_.direction(), // z-axis
        rotorGeometry_.psiRef()     // x-axis
    );

    cylCS_ = coordSystem::cylindrical // local Cartesian to cylindrical
        (
            rotorGeometry_.center(),    // centerd to local
            rotorGeometry_.direction(), // z-axis
            rotorGeometry_.psiRef()     // x-axis
        );

    carToCylCS_ = coordSystem::cylindrical(
        vector(0, 0, 0), // same center
        vector(0, 0, 1), // same z-axis
        vector(1, 0, 0)  // same x-axis
    );
}
void rotorDiscrete::writeArea(word propName, const fvMesh &mesh) const
{
    volScalarField areainfo(
        IOobject(
            propName + ":diskArea",
            mesh.time().timeName(),
            mesh),
        mesh,
        dimensionedScalar(dimArea, Zero));
    forAll(rotorCells_, i)
    {
        areainfo[rotorCells_[i].celli()] = rotorCells_[i].area();
    }

    areainfo.write();
}
tensor rotorDiscrete::bladeLocalFromPoint(const point &localPoint) const
{
    // z- up, y -outwards from center, x perpendicular y,z (leading edge to trailing edge)
    point global, origin;

    global = cylCS_.globalPosition(localPoint);
    origin = cylCS_.origin();

    tensor rotTensor(cylCS_.R());

    // z-axis
    const vector ax3 = rotTensor.col<2>(); // == e3 (already normalized)

    // y-axis (radial direction)
    vector ax2(global - origin);

    ax2.removeCollinear(ax3);

    const scalar magAxis2(mag(ax2));

    // Trap zero size and colinearity
    if (magAxis2 < SMALL)
    {
        return rotTensor;
    }

    ax2 /= magAxis2; // normalise

    // Replace with updated local axes

    rotTensor.col<0>(ax2 ^ ax3);
    rotTensor.col<1>(ax2);

    return rotTensor;
}

void rotorDiscrete::createMeshIntersect
(
    const rotorFvMeshSel& rotorFvMeshSel, 
    List<point> &vertex, 
    List<label>& cellidx, 
    List<List<label>> &cellPoints
)
{
    // Selected rotor radius (real used, no from mesh)
    const scalar radius = rotorGeometry_.radius();
    const scalar sqrRadius = radius * radius;
    const scalar idealArea = constant::mathematical::pi * sqrRadius;

    // List ref.
    const vectorField &cellCenter = rotorFvMeshSel.mesh().C();
    const labelListList &cellEdges = rotorFvMeshSel.mesh().cellEdges();
    const edgeList &edges = rotorFvMeshSel.mesh().edges();
    const vectorField &points = rotorFvMeshSel.mesh().points();

    // Reset
    vertex.resize(0);

    cellPoints.resize(cellidx.size());
    List<label> validCells(cellidx.size());
    label nValid = 0;
    vector n = rotorGeometry_.direction();
    point p0 = rotorGeometry_.center();

    forAll(cellidx, i)
    {
        label celli = cellidx[i];                   // Get cell idx
        labelList edgesInCell = cellEdges[celli]; // Get cell edges
        forAll(edgesInCell, j)
        {
            label edgei = edgesInCell[j];
            edge e = edges[edgei];
            // Line points
            point l0 = points[e.a()];
            point l1 = points[e.b()];

            vector l = l1 - l0;

            scalar den = l.inner(n);
            if (den == 0)
            {
                // no intersection
                // Or contained
            }
            else
            {
                scalar d = ((p0 - l0).inner(n)) / den;
                if (d >= -0.0 && d <= 1.0) // between l0 and l1
                {
                    vector intersect = l0 + l * d;
                    vector local =carCS_.localPosition(intersect);
                    cellPoints[nValid].append(vertex.size());
                    local.z()=0;
                    vertex.append(local);
                }
            }
        }
        if(cellPoints[nValid].size()>=3) //valid cut
        {
            validCells[nValid]=celli;
            nValid++;
        }
        else
        {
            cellPoints[nValid].clear();
        }
    }
    cellPoints.resize(nValid);
    validCells.resize(nValid);
    cellidx=validCells;
}

void rotorDiscrete::selectInnerCells(const rotorFvMeshSel &rotorFvMeshSel, List<label> &cellidx)
{
    List<label> cellSel = rotorFvMeshSel.cells();
    // Selected rotor radius (real used, no from mesh)
    const scalar radius = rotorGeometry_.radius();
    const scalar radiusSqr = radius * radius;

    const labelListList& cellVertex = rotorFvMeshSel.mesh().cellPoints(); 
    const pointField& meshPoints = rotorFvMeshSel.mesh().points();  

    cellidx.resize(cellSel.size());
    label nCell = 0;
    auto isInsideRotor = delaunayTriangulation::circularRegion(radius);
    
    forAll(cellSel,i)
    {
        label celli = cellSel[i];
        if(!includeVertex_)
        {
            //Include cell if the centroid proyects inside the rotor
            vector cellCentroid = rotorFvMeshSel.mesh().C()[celli];
            vector localPos = carCS_.localPosition(cellCentroid);
            localPos.z()=0; //set on rotor plane 

            if( isInsideRotor(localPos))
            {
                cellidx[nCell]=celli;
                nCell++;
            }
        }
        else
        {
            //Include if any cell vertex lies inside the rotor
            forAll(cellVertex[celli],j)
            {
                vector point = meshPoints[cellVertex[celli][j]];
                vector localPos = carCS_.localPosition(point);
                localPos.z()=0; //set on rotor plane 

                if(isInsideRotor(localPos))
                {
                    cellidx[nCell]=celli;
                    nCell++;
                    break;
                }
            }
        }
    }
    cellidx.resize(nCell); //nCell <= cellSel.size()
}

bool rotorDiscrete::cutWithCircle
    (
        List<point> &vertex,
        List<label> &cell,
        const std::function<bool(point)> &isInRegion,
        const std::function<void(point, point, point &, bool)> &findIntersection
    )
{
    // Create temp poligon
    List<label> poli = cell;
    // Add points inside domain or intersection with domain
    List<label> insidePoli;

    // List of visited points
    List<bool> visited(poli.size(), false);

    bool anyInside = false;
    forAll(poli, w)
    {
        if (isInRegion(vertex[poli[w]]))
        {
            anyInside = true;
        }
    }
    if (!anyInside)
    {
        cell.clear();
        return false;
    }
    for (label j = 0; j < poli.size(); j++)
    {

        point po;
        po = vertex[poli[j]];
        if (!isInRegion(po))
        {
            if (!visited[j])
            {
                point p1, p2, int1, int2, pout1, pout2;
                label i1, i2, o1, o2;

                List<label> to_remove;

                i1 = j;
                i2 = j;
                do
                {
                    visited[i1] = true;
                    o1 = i1;
                    i1 = i1 == 0 ? poli.size() - 1 : i1 - 1;
                    p1 = vertex[poli[i1]];
                } while (!isInRegion(p1));

                // Get next point indexo1 = (i1+1)%poli.size();

                do
                {
                    visited[i2] = true;
                    o2 = i2;
                    i2 = (i2 + 1) % poli.size();
                    p2 = vertex[poli[i2]];

                } while (!isInRegion(p2));

                if (o1 == o2)
                {
                    findIntersection(vertex[poli[o1]], p1, int1, true);
                    findIntersection(vertex[poli[o1]], p2, int2, true);

                    delaunayTriangulation::refineBorder(vertex, insidePoli, int1, int2, refinementLevel_, isInRegion, findIntersection);
                }
                else
                {
                    pout1 = vertex[poli[o1]];
                    pout2 = vertex[poli[o2]];

                    findIntersection(pout1, p1, int1, true);
                    findIntersection(pout2, p2, int2, true);

                    delaunayTriangulation::refineBorder(vertex, insidePoli, int1, int2, refinementLevel_, isInRegion, findIntersection);
                }
            }
        }
        else
        {
            insidePoli.append(poli[j]);
        }
    }
    cell = insidePoli;
    return true;
}

void rotorDiscrete::createFromData(const List<point> &vertex, const List<point> &centers,const List<label>& cellidx, List<List<label>> &cells)
{
    // Add vertex and center points to local data array
    carPoints_ = vertex;
    //carPoints_.resize(vertex.size() + centers.size());

    // Integration points enable flag
    integrationPoints_.resize(carPoints_.size(), false);

    //rotorCells_.resize(centers.size());

    // Selected rotor radius (real used, no from mesh)
    const scalar radius = rotorGeometry_.radius();
    const scalar sqrRadius = radius * radius;
    const scalar idealArea = constant::mathematical::pi * sqrRadius;


    auto isInsideRotor = delaunayTriangulation::circularRegion(radius);
    auto intersectRotor = delaunayTriangulation::intersectCircle(radius);
    area_ = 0.0;
    scalar ta = 0.0;
    rotorCells_.resize(cellidx.size());
    label goodIdx = 0;
    // add cell centers to rotor points and create cells
    forAll(cellidx, i)
    {
        label celli = cellidx[i]; // Now indexing is from cellis (used cells)
        delaunayTriangulation::sortCounterClockwise(carPoints_,cells[i]);
        if(!cutWithCircle(carPoints_,cells[i],isInsideRotor,intersectRotor))
        {
            continue;
        }
        

        vector c{0, 0, 0};
        // Centers the cell_center in the middle of the cell
        if (correctCenters_)
        {
            delaunayTriangulation::correctCenter(carPoints_,cells[i],c);
        }
        // Keeps the original cell_center
        else
        {
            //Correct only the centers outside the cell
            if(!delaunayTriangulation::isInsideCell(carPoints_,cells[i],centers[i]))
            {
                Warning<<"Center "<<i<<" :"<<centers[i] <<" lays outside the cell"<<endl;
                Warning<<"Size: "<<cells[i].size()<<endl;
                
                delaunayTriangulation::correctCenter(carPoints_,cells[i],c);
                if(!delaunayTriangulation::isInsideCell(carPoints_,cells[i],c))
                { 
                    if(!delaunayTriangulation::isConvex(carPoints_,cells[i]))
                    {
                        Warning<<"Poligon is not convex"<<endl;
                    }
                    Info<<"Point: "<<c<<endl<<endl;

                    FatalErrorInFunction
                        <<"Cell nº "<<i<<"is not valid and could not fix"
                        <<endl;
                }
            }
            else
            {
                c = centers[i];
            }
        }

        carPoints_.append(c); // add center
        label centerIdx = carPoints_.size()-1;
        
        rotorCells_.set(goodIdx,rotorCell::New(integrationMode_, centerIdx, cells[i], carPoints_, celli));

        integrationPoints_.resize(carPoints_.size(), false);
        rotorCells_[goodIdx].updateIntegrationList(integrationPoints_);

        area_ += rotorCells_[goodIdx].area();
        goodIdx++;
    }

    rotorCells_.resize(goodIdx);

    totalArea_ = area_;
    reduce(totalArea_, sumOp<scalar>());

    Info << "Total disk area: " << totalArea_ << endl;
    Info << "Ideal disk area: " << idealArea << endl;
    Info << "Area error = " << (idealArea - totalArea_) / idealArea * 100.0 << "%" << endl;

    // TODO: ?Add aditional integration points

    // Add cyl coord system and local blade tensor

    cylPoints_.resize(carPoints_.size());
    localBlade_.resize(carPoints_.size());
    forAll(carPoints_, i)
    {
        cylPoints_[i] = carToCylCS_.localPosition(carPoints_[i]);
        localBlade_[i] = this->bladeLocalFromPoint(cylPoints_[i]);
    }
}
void rotorDiscrete::writePythonPlotter(word process)
{
    std::string x_string = "x = [";
    std::string y_string = "y = [";
    std::string xc_string = "xc = [";
    std::string yc_string = "yc = [";

    for (label i = 0; i < rotorCells_.size(); i++)
    {

        vector center = rotorCells_[i].centerPosition();
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

    for (label i = 0; i < carPoints_.size(); i++)
    {
        x_string += std::to_string(carPoints_[i].x());

        y_string += std::to_string(carPoints_[i].y());

        if (i != carPoints_.size() - 1)
        {
            x_string += ",";
            y_string += ",";
        }
    }

    x_string += "]";
    y_string += "]";

    std::string tri_str = "tri = [";
    for (label i = 0; i < rotorCells_.size(); i++)
    {
        const auto& vor = rotorCells_[i].vertex();
        if (vor.size() == 0)
            continue;
        tri_str += "[";
        for (label j = 0; j < vor.size(); j++)
        {
            tri_str += std::to_string(vor[j]);
            tri_str += ",";
        }
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
    pyplot += "for i in range(len(tri)):\n";
    pyplot += "\txp = np.zeros(len(tri[i]))\n";
    pyplot += "\typ = np.zeros(len(tri[i]))\n";
    pyplot += "\tfor j in range(len(tri[i])):\n";
    pyplot += "\t\txp[j]=x[tri[i][j]]\n";
    pyplot += "\t\typ[j]=y[tri[i][j]]\n";
    // pyplot+="\tif(len(tri[i])>2):\n";
    pyplot += "\tplot.plot(xp,yp)\n";
    pyplot += "\n";
    pyplot += "plot.plot(xc,yc,marker='o',linewidth = 0)\n";
    pyplot += "plot.show()\n";
    std::ofstream file("triangulation" + std::to_string(Pstream::myProcNo()) + ".py", std::ios::out);
    file << pyplot;
    file.close();
}
void rotorDiscrete::fromMesh(const rotorFvMeshSel &rotorFvMeshSel, word integration, bool correctCenters, label refinementLevel,word discreteMethod)
{
    Info << "Building rotor Discrete from mesh" << endl;
    discreteMode_ = discreteMode::dmMesh;
    integrationMode_ = integration;
    refinementLevel_ = refinementLevel;
    discreteMethod_=discreteMethod;
    correctCenters_=correctCenters;

    // Selected rotor radius (real used, no from mesh)
    const scalar radius = rotorGeometry_.radius();
    const scalar sqrRadius = radius * radius;
    const scalar idealArea = constant::mathematical::pi * sqrRadius;

    //
    // List ref.
    const labelList &cells = rotorFvMeshSel.cells();
    const vectorField &cellCenter = rotorFvMeshSel.mesh().C();

    List<point> centers;
    List<label> cellis;

    Info<<"Mesh selection: "<<rotorFvMeshSel.cells().size()<<endl;
    selectInnerCells(rotorFvMeshSel,cellis);
    Info<<"Selected inner cells: "<<cellis.size()<<endl;

    List<List<label>> cellPoints;
    List<point> vertex;

    if(discreteMethod=="voronoid")
    {
        //Get cell centers
        centers.resize(cellis.size());
        forAll(centers,i)
        {
            centers[i]=carCS_.localPosition(cellCenter[cellis[i]]);
            centers[i].z()=0;
        }
        this->createMeshVoronoid(vertex,centers,cellPoints);
    }
    else
    {
        this->createMeshIntersect(rotorFvMeshSel,vertex,cellis,cellPoints);

        Info<<"After intersection: "<<cellis.size()<<endl;
        centers.resize(cellis.size());
        forAll(centers,i)
        {
            centers[i]=carCS_.localPosition(cellCenter[cellis[i]]);
            centers[i].z()=0;
        }
        
    }

    this->createFromData(vertex,centers,cellis,cellPoints);

    this->writePythonPlotter();
}


void rotorDiscrete::createMeshVoronoid
    (
        List<point> &vertex, 
        List<point> &centers,
        List<List<label>> &cells
    )
{

    // TODO
    // integration, correct centers, and refinement level shouldnt be provided to this fuction, should be in the dictionary!!!
    Info << "Creating voronoid discretization" << endl;

    // Gather ncell information
    List<label> ncellis(Pstream::nProcs(), 0);
    ncellis[Pstream::myProcNo()] = centers.size();
    reduce(ncellis, sumOp<labelList>());
    label totalCells = sum(ncellis);


    List<label> coreIdx(ncellis.size());
    label sum = 0;
    forAll(coreIdx, i)
    {
        coreIdx[i] = sum;
        sum += ncellis[i];
    }

    List<point> allCenter(totalCells, vector(0, 0, 0));

    auto it = centers.begin();
    auto dstIt = allCenter.begin() + coreIdx[Pstream::myProcNo()];
    while (it != centers.end())
    {
        (*dstIt) = (*it);
        it++;
        dstIt++;
    }
    reduce(allCenter, sumOp<vectorList>());

    List<List<label>> voroCells;

    scalar radius = rotorGeometry_.radius();

    // Create voronoid diagram of proyected centers for 2D meshing
    delaunayTriangulation::Voronoid(
        allCenter,
        vertex,
        voroCells,
        refinementLevel_,
        delaunayTriangulation::circularRegion(radius),
        delaunayTriangulation::intersectCircle(radius)
    );

    // Add vertex points
    cells.resize(centers.size());

    forAll(cells, i)
    {
        label globalIdx = i + coreIdx[Pstream::myProcNo()];
        cells[i]=voroCells[globalIdx];
    }

}
bool rotorDiscrete::read(const dictionary &dict)
{
    return false;
}
}