#include "meshGrid.H"
#include "cellSet.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "unitConversion.H"
#include "rotorGeometry.H"
#include "geometry.H"
#include <fstream>
#include "syncTools.H"
#include "meshCell.H"
#include "GatherList.H"
namespace Foam
{
defineTypeNameAndDebug(meshGrid, 0);
addToRunTimeSelectionTable(rotorGrid,meshGrid,dictionary);

const Enum
<
    meshGrid::discreteMethod
>
meshGrid::discreteMethodNames_
({
    {meshGrid::discreteMethod::voronoid, "voronoi"},
    {meshGrid::discreteMethod::intersection, "intersection"},
    {meshGrid::discreteMethod::proyection, "proyection"}
});

meshGrid::meshGrid
(
    const dictionary& dict,
    const rotorGeometry& geometry,
    const rotorFvMeshSel &rotorFvMeshSel,
    const bladeModelS* bladeModel,
    scalar nBlades
)
 : rotorGrid(dict,geometry,rotorFvMeshSel,nBlades)
{
    this->read(dict);
    assignFvCells();
    build();
}


void meshGrid::assignFvCells()
{
    this->fromMesh();
}


void meshGrid::createProyection(const List<label>& cells)
{
    vector axis = rotorGeometry_.direction();
    
    const auto& mesh_ = meshSel_.mesh(); 
    scalarField area(cells.size(),0.0);
    cells_.resize(cells.size());

    static const scalar tol = 0.8;

    const label nInternalFaces = mesh_.nInternalFaces();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const vectorField& Sf = mesh_.Sf();
    const scalarField& magSf = mesh_.magSf();

    vector n = Zero;

    // Calculate cell addressing for selected cells
    labelList cellAddr(mesh_.nCells(), -1);
    labelUIndList(cellAddr, cells) = identity(cells.size());
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
                area[own] += magSf[facei];
                n += Sf[facei];
            }
        }
        else if ((own == -1) && (nbr != -1))
        {
            vector nf = Sf[facei]/magSf[facei];

            if ((-nf & axis) > tol)
            {
                area[nbr] += magSf[facei];
                n -= Sf[facei];
            }
        }
    }

    // Add boundary contributions
    forAll(pbm, patchi)
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
                    area[own] += magSfp[j];
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
                    area[own] += magSfp[j];
                    n += Sfp[j];
                }
            }
        }
    }

    //Create cells
    label goodIdx = 0;
    const vectorField& meshCentroid = mesh_.C();


    //keep only the valid cells
    List<label> goodCell(area.size());
    List<scalar> goodArea(area.size());
    List<vector> goodCenter(area.size());
    
    const scalar radius = rotorGeometry_.radius();
    const scalar minradius = rotorGeometry_.innerRadius();
    auto isInsideCorona = util::geometry::coronaRegion(radius,minradius);


    forAll(area,i)
    {
        vector cellCenter = meshCentroid[cells[i]];
        cellCenter = rotorGeometry_.cylindricalCS().localPosition(cellCenter);
        cellCenter.z()=0;
        cellCenter = rotorGeometry_.cartesianToCylindrical().globalPosition(cellCenter);

        if(area[i]>SMALL && isInsideCorona(cellCenter))
        {
            goodCenter[goodIdx] = cellCenter;
            
            goodCell[goodIdx]=cells[i];
            goodArea[goodIdx]=area[i];
            ++goodIdx;
        }
    }
    goodCell.resize(goodIdx);
    goodArea.resize(goodIdx);
    goodCenter.resize(goodIdx);


    //Gather par Data
    label totalCell = goodIdx;
    reduce(totalCell,sumOp<label>());

    List<label> cellSize(Pstream::nProcs(),0);
    cellSize[Pstream::myProcNo()]=goodIdx;
    reduce(cellSize,sumOp<labelList>());

    List<List<label>> parCells(Pstream::nProcs());
    List<List<scalar>> parArea(Pstream::nProcs());
    List<List<vector>> parCenter(Pstream::nProcs());

    forAll(parCells,i)
    {
        parCells[i].resize(cellSize[i],0);
        parArea[i].resize(cellSize[i],0);
        parCenter[i].resize(cellSize[i],Zero);
    }
    parCells[Pstream::myProcNo()]=goodCell;
    parArea[Pstream::myProcNo()]=goodArea;
    parCenter[Pstream::myProcNo()]=goodCenter;
   
    forAll(parCells,i)
    {
        reduce(parCells[i],sumOp<labelList>());
        reduce(parArea[i],sumOp<scalarList>());
        reduce(parCenter[i],sumOp<vectorList>());
    }

    cells_.resize(totalCell);
    label cellCount = 0;
    forAll(parCells,i)
    {
        auto& localCell = parCells[i];
        auto& localArea = parArea[i];
        auto& localCenter = parCenter[i];
        forAll(localCell,j)
        {
            label cellidx = localCell[j];
            scalar cellArea = localArea[j];
            if(Pstream::myProcNo()!=i)
            {
                //No index in this core
                cellidx=-1;
            }
            area_+=cellArea;
            meshCell* newcell = new meshCell(rotorGeometry_,cellidx,nBlades_,cellArea,localCenter[j]);
            cells_.set(cellCount,newcell);
            ++cellCount;
        }
    }
    cells_.resize(cellCount);

}

void meshGrid::createMeshIntersect
(
    List<point> &vertex, 
    List<label>& cellidx, 
    List<List<label>> &cellPoints,
    List<bool> &isThisCore
)
{
    // List ref.
    const labelListList &cellEdges = meshSel_.mesh().cellEdges();
    const edgeList &edges = meshSel_.mesh().edges();
    const vectorField &points = meshSel_.mesh().points();

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
                    vector local =rotorGeometry_.cartesianCS().localPosition(intersect);
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

    List<label> vertexIdx;
    vertex = util::GatherList(vertex,vertexIdx);

    forAll(cellPoints,i)
    {
        forAll(cellPoints[i],j)
        {
            cellPoints[i][j]+=vertexIdx[Pstream::myProcNo()];
        }
    }

    List<label> cellPointsIdx;
    cellPoints = util::GatherListList<label>(cellPoints,cellPointsIdx);

    isThisCore.resize(cellPoints.size(),false);
    forAll(validCells,i)
    {
        label idx = i + cellPointsIdx[Pstream::myProcNo()];
        isThisCore[idx]=true;
    }
}

void meshGrid::selectInnerCells(List<label> &cellidx)
{
    List<label> cellSel = meshSel_.cells();
    // Selected rotor radius (real used, no from mesh)
    const scalar radius = rotorGeometry_.radius();

    const labelListList& cellVertex = meshSel_.mesh().cellPoints(); 
    const pointField& meshPoints = meshSel_.mesh().points();  

    cellidx.resize(cellSel.size());
    label nCell = 0;
    auto isInsideRotor = util::geometry::circularRegion(radius);
    
    forAll(cellSel,i)
    {
        label celli = cellSel[i];
        if(!includeVertex_)
        {
            //Include cell if the centroid proyects inside the rotor
            vector cellCentroid = meshSel_.mesh().C()[celli];
            vector localPos = rotorGeometry_.cylindricalCS().localPosition(cellCentroid);
            localPos.z()=0; //set on rotor plane 
            localPos = rotorGeometry_.cartesianToCylindrical().globalPosition(localPos);
            if(isInsideRotor(localPos))
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
                vector localPos = rotorGeometry_.cylindricalCS().localPosition(point);
                localPos.z()=0; //set on rotor plane 
                localPos = rotorGeometry_.cartesianToCylindrical().globalPosition(localPos);

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

bool meshGrid::cutWithCircle
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

                    util::geometry::refineBorder(vertex, insidePoli, int1, int2, refinementLevel_, isInRegion, findIntersection);
                }
                else
                {
                    pout1 = vertex[poli[o1]];
                    pout2 = vertex[poli[o2]];

                    findIntersection(pout1, p1, int1, true);
                    findIntersection(pout2, p2, int2, true);

                    util::geometry::refineBorder(vertex, insidePoli, int1, int2, refinementLevel_, isInRegion, findIntersection);
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

void meshGrid::createFromData
(
    List<point> &vertex,
    const List<point> &centers,
    const List<label>& cellidx,
    List<List<label>> &cells,
    List<bool>& isThisCore
)
{    
    // Selected rotor radius (real used, no from mesh)
    const scalar radius = rotorGeometry_.radius();
    const scalar minradius = rotorGeometry_.innerRadius();
    label goodIdx = 0;
    label coreCount = 0;
    label layOut = 0;

    auto isInsideRotor = util::geometry::circularRegion(radius);
    auto intersectRotor = util::geometry::intersectCircle(radius);
    auto isInsideCorona = util::geometry::coronaRegion(radius,minradius);
    cells_.resize(cells.size());
    // add cell centers to rotor points and create cells
    forAll(cells, i)
    {

        util::geometry::sortCounterClockwise(vertex,cells[i]);
        if(!cutWithCircle(vertex,cells[i],isInsideRotor,intersectRotor))
        {
            continue;
        }
    
        vector c{0, 0, 0};
        // Centers the cell_center in the middle of the cell
        if (correctCenters_)
        {
            util::geometry::correctCenter(vertex,cells[i],c);
        }
        // Keeps the original cell_center
        else
        {
            //Correct only the centers outside the cell
            if(!util::geometry::isInsideCell(vertex,cells[i],centers[i]))
            {
                layOut++;
                
                util::geometry::correctCenter(vertex,cells[i],c);
                if(!util::geometry::isInsideCell(vertex,cells[i],c))
                { 
                    if(!util::geometry::isConvex(vertex,cells[i]))
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


        label celli;
        if(isThisCore[i])
        {
            celli = cellidx[coreCount]; // Now indexing is from cellis (used cells)
            ++coreCount;
        }
        else
        {
            celli=-1;
        }
        if(!isInsideCorona(c))
        {
            continue;
        }
        meshCell* newcell = new meshCell(rotorGeometry_,celli, nBlades_,vertex,cells[i],&c);
        cells_.set(goodIdx,newcell);
        area_+=newcell->area();
        goodIdx++;
    }
    cells_.resize(goodIdx);

    if(layOut>0)
    {
        Warning
            <<layOut<<" cells centers lay outsid rotorCell and have been corrected"
            <<endl;
    }
}
void meshGrid::writePythonPlotter(word outputName)
{

    if(discreteMethod_==discreteMethod::proyection)
    {
        //Not enabled for this one
        return;
    }
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
        meshCell* cell = dynamic_cast<meshCell*>(cells_.get(i));

        auto& points = cell->localCartesianPoints();
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
void meshGrid::fromMesh()
{
    Info<<endl;
    Info << "Building rotor Discrete from mesh:" << endl;
    Info.stream().incrIndent();

    // Selected rotor radius (real used, no from mesh)
    const scalar radius = rotorGeometry_.radius();
    const scalar minradius = rotorGeometry_.innerRadius();

    const scalar sqrRadius = radius * radius;
    const scalar sqrMinradius = minradius * minradius;

    const scalar idealArea = constant::mathematical::pi * (sqrRadius-sqrMinradius);

    List<point> centers;
    List<label> cellis;

    selectInnerCells(cellis);

    if(discreteMethod_ == discreteMethod::proyection)
    {
        createProyection(cellis);
    }
    else
    {
        const vectorField &cellCenter = meshSel_.mesh().C();

        List<point> vertex;
        List<List<label>> cellPoints;
        List<bool> isThisCore;

        if(discreteMethod_==discreteMethod::voronoid)
        {
            //Get cell centers
            centers.resize(cellis.size());
            forAll(centers,i)
            {
                centers[i]=rotorGeometry_.cylindricalCS().localPosition(cellCenter[cellis[i]]);
                centers[i].z()=0;
                centers[i]=rotorGeometry_.cartesianToCylindrical().globalPosition(centers[i]);
            }
            this->createMeshVoronoid(vertex,centers,cellPoints,isThisCore);
        }
        else
        {
            this->createMeshIntersect(vertex,cellis,cellPoints,isThisCore);

            centers.resize(cellis.size());
            forAll(centers,i)
            {
                centers[i]=rotorGeometry_.cylindricalCS().localPosition(cellCenter[cellis[i]]);
                centers[i].z()=0;
                centers[i]=rotorGeometry_.cartesianToCylindrical().globalPosition(centers[i]);
            }
            labelList centIdx;
            centers = util::GatherList<vector>(centers,centIdx);
            
        }

        this->createFromData(vertex,centers,cellis,cellPoints,isThisCore);
    }


    //Compute area
    indent(Info) << "- Total disk area: " << area_ << endl;
    indent(Info) << "- Area error = " << (idealArea - area_) / idealArea * 100.0 << "%" << endl;

    Info.stream().decrIndent();

    this->writeArea("propeller1:");
}


void meshGrid::createMeshVoronoid
    (
        List<point> &vertex, 
        List<point> &centers,
        List<List<label>> &cells,
        List<bool>& isThisCore
    )
{

    // Gather ncell information
    List<label> coreIdx;
    List<vector> allCenter = util::GatherList<vector>(centers,coreIdx);

    scalar radius = rotorGeometry_.radius();
    // Create voronoid diagram of proyected centers for 2D meshing
    util::geometry::Voronoid(
        allCenter,
        vertex,
        cells,
        refinementLevel_,
        util::geometry::circularRegion(radius),
        util::geometry::intersectCircle(radius)
    );

    // Add vertex points
    isThisCore.resize(0);
    isThisCore.resize(cells.size(),false);

    //Set the core cells
    forAll(centers, i)
    {
        label globalIdx = i + coreIdx[Pstream::myProcNo()];
        isThisCore[globalIdx]=true;
    }
    centers=allCenter;
}

bool meshGrid::read(const dictionary &dict)
{
    discreteMethod_ = discreteMethodNames_.getOrDefault("discreteMethod",dict,discreteMethod::voronoid);

    //For voronoid method vertex cells are not needed, but for intersection is recommended
    bool defaultVertex = discreteMethod_ == discreteMethod::intersection;
    includeVertex_ = dict.getOrDefault<bool>("includeIfVertex",defaultVertex);

    refinementLevel_ = dict.getOrDefault<label>("borderRefinement",0);

    correctCenters_ = dict.getOrDefault<bool>("correctCenters",false);
    
    Info<<endl;    
    Info << "Reading rotorGrid dict:" << endl;
    Info.stream().incrIndent();
    indent(Info)<<"- Discrete method: "<<discreteMethodNames_.get(discreteMethod_)<<endl;
    indent(Info)<<"- Border refinement: "<<refinementLevel_<<endl;
    indent(Info)<<"- Correct centers: "<<correctCenters_<<endl;
    indent(Info)<<"- Include if vertex inside: "<<includeVertex_<<endl;
    Info.stream().decrIndent();

    return true;
}


}
