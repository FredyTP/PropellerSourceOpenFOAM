#include "rotorDiscrete.H"
#include "cellSet.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "unitConversion.H"
#include "rotorGeometry.H"
#include "delaunayTriangulation.H"
#include <fstream>
#include "rotorTriCell.H"
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
void rotorDiscrete::createMeshOld(const rotorFvMeshSel& rotorFvMeshSel)
{
    vector axis = rotorGeometry_.direction();

    const auto& mesh_ = rotorFvMeshSel.mesh(); 
    const auto& cells_ = rotorFvMeshSel.cells();
    scalarField area_(cells_.size(),1.0);

    static const scalar tol = 0.8;

    const label nInternalFaces = mesh_.nInternalFaces();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const vectorField& Sf = mesh_.Sf();
    const scalarField& magSf = mesh_.magSf();

    vector n = Zero;

    // Calculate cell addressing for selected cells
    labelList cellAddr(mesh_.nCells(), -1);
    labelUIndList(cellAddr, cells_) = identity(cells_.size());
    labelList nbrFaceCellAddr(mesh_.nBoundaryFaces(), -1);
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label facei = pp.start() + i;
                label nbrFacei = facei - nInternalFaces;
                label own = mesh_.faceOwner()[facei];
                nbrFaceCellAddr[nbrFacei] = cellAddr[own];
            }
        }
    }

    // Correct for parallel running
    syncTools::swapBoundaryFaceList(mesh_, nbrFaceCellAddr);

    //Only computes de area of the cellSet boundary which has the same normal than rotor normal, because only takes
    //Into account the face that have no owner or no neighbour in the selection

    // Add internal field contributions
    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        const label own = cellAddr[mesh_.faceOwner()[facei]];
        const label nbr = cellAddr[mesh_.faceNeighbour()[facei]];

        if ((own != -1) && (nbr == -1))
        {
            vector nf = Sf[facei]/magSf[facei];

            if ((nf & axis) > tol)
            {
                area_[own] += magSf[facei];
                n += Sf[facei];
            }
        }
        else if ((own == -1) && (nbr != -1))
        {
            vector nf = Sf[facei]/magSf[facei];

            if ((-nf & axis) > tol)
            {
                area_[nbr] += magSf[facei];
                n -= Sf[facei];
            }
        }
    }


    // Add boundary contributions
    /*forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];
        const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];

        if (pp.coupled())
        {
            forAll(pp, j)
            {
                const label facei = pp.start() + j;
                const label own = cellAddr[mesh_.faceOwner()[facei]];
                const label nbr = nbrFaceCellAddr[facei - nInternalFaces];
                const vector nf = Sfp[j]/magSfp[j];

                if ((own != -1) && (nbr == -1) && ((nf & axis) > tol))
                {
                    area_[own] += magSfp[j];
                    n += Sfp[j];
                }
            }
        }
        else
        {
            forAll(pp, j)
            {
                const label facei = pp.start() + j;
                const label own = cellAddr[mesh_.faceOwner()[facei]];
                const vector nf = Sfp[j]/magSfp[j];

                if ((own != -1) && ((nf & axis) > tol))
                {
                    area_[own] += magSfp[j];
                    n += Sfp[j];
                }
            }
        }
    }*/

    volScalarField areaIO
    (
        IOobject
        (
            "oldArea",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimArea, Zero)
    );
    UIndirectList<scalar>(areaIO.primitiveField(), cells_) = area_;

    Info<< " writing field " << areaIO.name()
        << endl;

    areaIO.write();
    
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
    // List ref.
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

    auto isInsideRotor = delaunayTriangulation::circularRegion(radius);
    auto intersectRotor = delaunayTriangulation::intersectCircle(radius);
    area_ = 0.0;
    rotorCells_.resize(cellidx.size());
    label goodIdx = 0;
    labelList layOut;
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
                layOut.append(cells[i]);
                
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
    // TODO: ?Add aditional integration points
    if(layOut.size()>0)
    {
        Warning
            <<layOut.size()<<" cells centers lay outsid rotorCell and have been corrected"
            <<endl;
    }
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
void rotorDiscrete::fromMesh(const rotorFvMeshSel &rotorFvMeshSel)
{
    Info<<endl;
    Info << "Building rotor Discrete from mesh:" << endl;
    Info.stream().incrIndent();

    discreteMode_ = discreteMode::dmMesh;

    // Selected rotor radius (real used, no from mesh)
    const scalar radius = rotorGeometry_.radius();
    const scalar sqrRadius = radius * radius;
    const scalar idealArea = constant::mathematical::pi * sqrRadius;

    this->createMeshOld(rotorFvMeshSel);
    //
    // List ref.
    const vectorField &cellCenter = rotorFvMeshSel.mesh().C();
    const scalarField& cellVol = rotorFvMeshSel.mesh().V();
    List<point> centers;
    List<label> cellis;

    selectInnerCells(rotorFvMeshSel,cellis);

    List<List<label>> cellPoints;
    List<point> vertex;

    if(discreteMethod_==discreteMethod::dmeVoronoid)
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

        centers.resize(cellis.size());
        forAll(centers,i)
        {
            centers[i]=carCS_.localPosition(cellCenter[cellis[i]]);
            centers[i].z()=0;
        }
        
    }

    this->createFromData(vertex,centers,cellis,cellPoints);

    //Join all core cell number
    label allCoreCells = rotorCells_.size();
    reduce(allCoreCells,sumOp<label>());

    //Obtain number of cells for each core
    labelList parNcells(Pstream::nProcs(),0);
    parNcells[Pstream::myProcNo()]=rotorCells_.size();

    reduce(parNcells,sumOp<labelList>());
     //Out total number of cells
    indent(Info)<<"- Total created rotorCells: "<<allCoreCells<<endl;
    indent(Info)<<"    - In each core: "<<parNcells<<endl;

    scalarList eachArea(Pstream::nProcs(),0);
    eachArea[Pstream::myProcNo()]=area_;
    reduce(eachArea,sumOp<scalarList>());

    indent(Info) << "- Total disk area: " << totalArea_ << endl;
    indent(Info) << "    - In each core: " << eachArea << endl;
    indent(Info) << "- Ideal disk area: " << idealArea << endl;
    indent(Info) << "- Area error = " << (idealArea - totalArea_) / idealArea * 100.0 << "%" << endl;

    Info.stream().decrIndent();

    this->writePythonPlotter();

    scalar maxRadius = geometry().radius();
    grid = rotorGrid(nRadial,nAzimutal,0.1*maxRadius,maxRadius);

    forAll(rotorCells_,i)
    {
        vector p = rotorCells_[i].centerPosition();
        vector polar = carToCylCS_.localPosition(p);
        label ir=0;
        label unused;
        label it=0;
        int result = regularInterpolation<scalar,scalar,1>::FindIndex(polar.x(),grid.radius(),ir,unused);
        regularInterpolation<scalar,scalar,1>::FindIndex(polar.y(),grid.theta(),it,unused);

        if(result == 1)
        {
            grid.cell(ir,it).addCelli(rotorCells_[i].celli(),cellVol[rotorCells_[i].celli()]);
        }
    }

    grid.build();
}


void rotorDiscrete::createMeshVoronoid
    (
        List<point> &vertex, 
        List<point> &centers,
        List<List<label>> &cells
    )
{

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

    integrationMode_ = rotorCell::integrationModeNames_.getOrDefault("integrationMode",dict,rotorCell::integrationMode::imCenter);
    discreteMethod_ = discreteMethodNames_.getOrDefault("discreteMethod",dict,discreteMethod::dmeVoronoid);

    //For voronoid method vertex cells are not needed, but for intersection is recommended
    bool defaultVertex = discreteMethod_ == discreteMethod::dmeIntersection;
    includeVertex_ = dict.getOrDefault<bool>("includeIfVertex",defaultVertex);

    refinementLevel_ = dict.getOrDefault<label>("borderRefinement",0);

    correctCenters_ = dict.getOrDefault<bool>("correctCenters",false);
    
    nRadial = dict.get<label>("nRadial");
    nAzimutal = dict.get<label>("nAzimutal");

    Info<<endl;    
    Info << "Reading rotor Discrete dict:" << endl;
    Info.stream().incrIndent();
    indent(Info)<<"- Discrete method: "<<rotorDiscrete::discreteMethodNames_.get(discreteMethod_)<<endl;
    indent(Info)<<"- Border refinement: "<<refinementLevel_<<endl;
    indent(Info)<<"- Correct centers: "<<correctCenters_<<endl;
    indent(Info)<<"- Integration mode: "<<rotorCell::integrationModeNames_.get(integrationMode_)<<endl;
    indent(Info)<<"- Include if vertex: "<<includeVertex_<<endl;
    Info.stream().decrIndent();

    return true;
}
}
