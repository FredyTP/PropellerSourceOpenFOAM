#include "rotorMesh.H"

#include "rotorMesh.H"
#include "cellSet.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "unitConversion.H"
#include "cellBitSet.H"
#include "cartesianCS.H"
#include "delaunayTriangulation.H"

namespace Foam
{

defineTypeNameAndDebug(rotorMesh,0);

const Enum
<
    rotorMesh::selectionMode
>
rotorMesh::selectionModeNames_
({
    {selectionMode::smNone, "none"},
    {selectionMode::smGeometry, "geometry"},
    {selectionMode::smCellSet, "cellSet"},
    {selectionMode::smCellZone, "cellZone"},
    {selectionMode::smGeometry2, "geometry2"}
});




rotorMesh::rotorMesh
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    selMode_(selectionMode::smNone),
    cells_(),
    area_(),
    diskArea_(NO_AREA),
    meshGeometry_(),
    cellsName_(),
    findClosestCenter_(false),
    correctGeometry_(false)
{

}

void rotorMesh::build(rotorGeometry& rotorGeometry)
{
    // Fix rotor center specification
    //Provide geometry data if not data specified or ask for correction
    this->clear();

    switch (selMode_)
    {
    case selectionMode::smGeometry :
    case selectionMode::smGeometry2 :
    {
        meshGeometry_= rotorGeometry; //Use as input data

        const bool rgReady = rotorGeometry.isRadiusReady() 
                && rotorGeometry.isCenterReady() 
                && rotorGeometry.isDirectionReady();

        //Check if geometry is enough to build mesh
        if(!rgReady)
        {
            FatalErrorInFunction
                << "Trying to select rotor mesh with incomplete geometry"
                <<exit(FatalError);
        }


        //Update to closest center if it's asked to (revisar esto)
        this->tryUpdateCenter(); 

        //Build local coordinate system
        localCartesianCS_ = coordSystem::cartesian
            (
                meshGeometry_.center(), //centerd to local
                meshGeometry_.direction(), //z-axis 
                meshGeometry_.psiRef()  //x-axis
            );

        if(selMode_==selectionMode::smGeometry)
        {
            this->createMeshSelection(); //Create cell selection from a disk
        }
        else
        {
            this->createMeshSelectionCilinder();
        }

        this->computeLocalPositions();

        Info<<"Proyected Area : "<<endl;
        this->computeCellsArea(); //

        Info<<"Voronoi Area: "<<endl;
        this->computeCellsAreaVoronoid();

        //Find "new" geometry --
        this->findRotorRadius();
        //Dont find rotor center because already selected cells from spec or aprox center
        //this->findRotorCenter();
        this->findRotorNormal(rotorGeometry);
        this->tryCorrectGeometry(rotorGeometry); //Update data if required

        if(meshGeometry_.radius()<rotorGeometry.radius())
        {
            FatalErrorInFunction
            <<"Mesh Selection radius is inferior to rotor radius"
            <<endl;
        }
    }
        break;
    case selectionMode::smCellZone :
    case selectionMode::smCellSet :

        //Check for input information
        if(!rotorGeometry.isDirectionReady())
        {
            FatalErrorInFunction
            <<"No defined disk normal for rotor mesh selection"
            <<exit(FatalError);
        }

        this->loadMeshSelection(); //Load mesh from cellSet

        //Find geometry---------
        this->findRotorCenter();
        this->findRotorNormal(rotorGeometry); //ROTOR NORMAL ALWAYS UPDATES ROTOR GEOMETRY NORMAL
        

        localCartesianCS_ = coordSystem::cartesian
            (
                meshGeometry_.center(), //centerd to local
                meshGeometry_.direction(), //z-axis 
                meshGeometry_.psiRef()  //x-axis
            );

        this->computeLocalPositions();
        this->findRotorRadius();

        Info<<"Proyected Area : "<<endl;
        this->computeCellsArea(); //

        Info<<"Voronoi Area: "<<endl;
        this->computeCellsAreaVoronoid();

        this->tryCorrectGeometry(rotorGeometry); //Use geometry as output
    break;

    default:
    break;
    }


    
    built_ = true;
}
void rotorMesh::clear()
{
    cells_.clear();
    area_.clear();
    diskArea_ = 0;
    built_ = false;
}
bool rotorMesh::read(const dictionary &dict)
{
    Info<<"Reading Rotor mesh config"<<endl;

    bool ok = true;

    selMode_ = selectionModeNames_.getOrDefault("selectionMode",dict,selectionMode::smGeometry);

    switch (selMode_)
    {
    case selectionMode::smGeometry :
    case selectionMode::smGeometry2 :
          
        findClosestCenter_ = dict.getOrDefault<bool>("closestCenter",false);
        Info<<"FindClosestCenter: "<< findClosestCenter_<<endl;

        //If rotorMesh is build from geometry, the resulting mesh geometry
        // may be used to update the provided geometry
        correctGeometry_ = dict.getOrDefault<bool>("correctGeometry",false);

        break;

    case selectionMode::smCellZone :
        ok &= dict.readEntry("cellZone",cellsName_);
        Info<< "Reading cellSet from: " << cellsName_<<endl;

        //If rotorMesh is from cellset geometry may not be provided
        correctGeometry_ = dict.getOrDefault<bool>("correctGeometry",false);

        break;
    case selectionMode::smCellSet :

        ok &= dict.readEntry("cellSet",cellsName_);
        Info<< "Reading cellSet from: " << cellsName_<<endl;

        //If rotorMesh is from cellset geometry may not be provided
        correctGeometry_ = dict.getOrDefault<bool>("correctGeometry",false);

        break;
    default:
        break;
    }



    return ok;
}
void rotorMesh::tryUpdateCenter()
{

    if(findClosestCenter_)
    {
        //Use closest cell centroid to center the plane used
        //to find rotor geometry, usually provides better results
        //in an oriented mesh
        vector oldCenter = meshGeometry_.center();
        label centercell = mesh_.findCell(oldCenter);
        if(centercell == -1)
        {
            //Maybe check before if point is inside of boundaries
            FatalErrorInFunction << "Rotor center is outside of mesh" 
            << exit(FatalError);
            
        }
        meshGeometry_.setCenter(mesh_.C()[centercell]);


        Info << "Using rotor center: " << meshGeometry_.center()
            << " instead of: " << oldCenter 
            << endl;
    }
}
void rotorMesh::tryCorrectGeometry(rotorGeometry& rotorGeometry)
{
    if(correctGeometry_)
    {
        rotorGeometry.setRadius(meshGeometry_.radius());
        rotorGeometry.setCenter(meshGeometry_.center());
        rotorGeometry.setDirection(meshGeometry_.direction());
        //Check for no refPsi direction (?)
    }
    else
    {
        //Add missing geometry information
        //Useful when selecting from mesh
        //and some parameters want to be specified like radius
        //when a mesh selection is bigger than the desired radius
        //And others want to be deduced like center or direction

        if(!rotorGeometry.isRadiusSet())
        {
            rotorGeometry.setRadius(meshGeometry_.radius());
        }
        if(!rotorGeometry.isCenterSet())
        {
            rotorGeometry.setCenter(meshGeometry_.center());
        }
        if(!rotorGeometry.isDirectionSet())
        {
            rotorGeometry.setCenter(meshGeometry_.direction());
        }
    }
}
void rotorMesh::computeLocalPositions()
{
    cellCenterLC_.resize(cells_.size());

    forAll(cells_,i)
    {
        label celli = cells_[i];

        vector cellCentroid = mesh_.C()[celli];
        vector localPos = localCartesianCS_.localPosition(cellCentroid);
        
        cellCenterLC_[i]= localPos;

    }
}
void rotorMesh::createMeshSelectionCilinder()
{
    const vector& rotorDir = meshGeometry_.direction();
    const vector& rotorCenter = meshGeometry_.center();
    const scalar& radius = meshGeometry_.radius();

    dictionary diskDict;
    diskDict.add("name","diskSelection");
    diskDict.add("type","cellSet");
    diskDict.add("action","add");
    diskDict.add("source","cylinderToCell");

    scalar width = 0.01;
    vector p1 = rotorCenter - width/2*rotorDir;
    vector p2 = rotorCenter + width/2*rotorDir;
    diskDict.add("p1",p1);
    diskDict.add("p2",p2);
    diskDict.add("radius",radius);

    dictionary selections;
    selections.add("disk",diskDict);

    bitSet selectedCells
    (
        cellBitSet::select(mesh_,selections,true)
    );
    cells_ = selectedCells.sortedToc();
}
void rotorMesh::createMeshSelection()
{

    const vector& rotorDir = meshGeometry_.direction();
    const vector& rotorCenter = meshGeometry_.center();
    const scalar& radius = meshGeometry_.radius();
    //TODO: Find better or improve algorithm to select rotor cells
    plane rotorPlane(rotorCenter,rotorDir);

    cuttingPlane cutPlane(rotorPlane,mesh_,false);

    labelList& planeCells = cutPlane.meshCells();

    //User squared radius to compute faster
    scalar radiusSqr = radius * radius;

    //Select only cells which centroid to the rotor center <= radius
    forAll(planeCells,i)
    {
        
        label celli = planeCells[i];

        vector cellCentroid = mesh_.C()[celli];
        vector localPos = localCartesianCS_.localPosition(cellCentroid);
        //localPos.z()=0; //set on rotor plane 
        scalar distanceSqr = magSqr(localPos);

        if( distanceSqr <= radiusSqr )
        {
            cells_.append(celli);
            cellCenterLC_.append(localPos);
        }
    }

    

    Info<<"Number of Selected cells: "<<cells_.size()<<endl;

}
void rotorMesh::loadMeshSelection()
{

    switch (selMode_)
    {    
    case selectionMode::smCellSet:
        //- Load from cell Set
        cells_ = cellSet(mesh_, cellsName_).sortedToc();

        break;
    case selectionMode::smCellZone:
    {
//- LOAD FROM CELLZONE
        Info<< indent
        << "- selecting cells using cellZones "
        << flatOutput(cellsName_) << nl;

        const auto& zones = mesh_.cellZones();

        label zoneId = zones.findIndex(cellsName_);

        if (zoneId == -1)
        {
            FatalErrorInFunction
                << "No matching cellZones: "
                << cellsName_ << nl
                << "Valid zones : "
                << flatOutput(zones.names()) << nl
                << "Valid groups: "
                << flatOutput(zones.groupNames())
                << nl
                << exit(FatalError);
        }
        //- Find more elegant way, searching for index twice
        cells_ = zones.selection(cellsName_,true).sortedToc();
        break;
    }
    case selectionMode::smGeometry :
    case selectionMode::smNone :
    default:
        break;
    }

    Info<<"Number of Selected cells: "<<cells_.size()<<endl;
    
}
void rotorMesh::computeCellsArea()
{
    const vector& rotorDir = meshGeometry_.direction();
    const scalar& radius = meshGeometry_.radius();

    scalar idealArea = constant::mathematical::pi * radius * radius;
    
    //set area size
    area_.resize(cells_.size());

    area_= 0.0;
    diskArea_ = 0.0;
    
    const label nInternalFaces = mesh_.nInternalFaces();
    const vectorField& Sf = mesh_.Sf();


    // Calculate cell addressing for selected cells
    labelList cellAddr(mesh_.nCells(), -1);
    labelUIndList(cellAddr, cells_) = identity(cells_.size());

    // Add internal field contributions
    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        const label own = cellAddr[mesh_.faceOwner()[facei]];
        const label nbr = cellAddr[mesh_.faceNeighbour()[facei]];

        if ((own != -1) && (nbr == -1))
        {
            area_[own] += 0.5*std::abs(Sf[facei] & rotorDir);
        }
        else if ((own == -1) && (nbr != -1))
        {
            area_[nbr] += 0.5*std::abs(Sf[facei] & rotorDir);
        }
    }

    forAll(area_,i)
    {
        diskArea_ += area_[i];
    }

    volScalarField areainfo
        (
            IOobject
            (
                "rotorDisk:area",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar(dimArea, Zero)
        );
    UIndirectList<scalar>(areainfo.primitiveField(), cells_) = area_;

    areainfo.write();
    Info<< "Ideal disk area: " << idealArea<<endl;
    Info<< "Disk Area: " << diskArea_ << endl;    

}

void rotorMesh::computeCellsAreaVoronoid()
{
    const vectorField& cellCenter = mesh_.C();
    const scalar radius = meshGeometry_.radius();
    scalar idealArea = constant::mathematical::pi * radius * radius;
    
    //set area size
    area_.resize(cells_.size());

    area_= 0.0;
    diskArea_ = 0.0;

    //Find cells centers
    List<point> centers(cells_.size());
    forAll(cells_,i)
    {
        label celli = cells_[i];
        centers[i] = localCartesianCS_.localPosition(cellCenter[celli]);
        centers[i].z() = 0; //proyect over plane
    }

    List<List<label>> voroCells;
    List<point> vertex;

    
    delaunayTriangulation::Voronoid
    (
        centers,
        vertex,    //new order no ref
        voroCells, // same order as sortedCenter
        delaunayTriangulation::circularRegion(radius),
        delaunayTriangulation::intersectCircle(radius)
    );

    // set vertex points
    rotorPoints_ = vertex; 

    rotorCells_.resize(cells_.size());
    rotorPoints_.resize(vertex.size() + centers.size());
    // add cell centers to rotor points and create cells
    forAll(voroCells,i)
    {
        label celli = cells_[i];
        rotorPoints_[i+vertex.size()] = centers[i];
        rotorCells_[i]=rotorCell(i+vertex.size(),voroCells[i],rotorPoints_,mesh_.V()[celli]);
    }

    //Find total area
    forAll(rotorCells_,i)
    {
        area_[i] = rotorCells_[i].area();
        diskArea_ += rotorCells_[i].area();
    }

    volScalarField areainfo
    (
        IOobject
        (
            "rotorDisk:area",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimArea, Zero)
    );
    UIndirectList<scalar>(areainfo.primitiveField(), cells_) = area_;

    areainfo.write();
    Info<< "Ideal disk area: " << idealArea<<endl;
    Info<< "Disk Area: " << diskArea_ << endl;    

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

void rotorMesh::findRotorCenter()
{
    const scalarField& cellVolume = mesh_.V();
    const vectorField& cellCenter = mesh_.C();

    vector newCenter(Zero);
    scalar volume = 0;
    
    forAll(cells_, i)
    {   
        label celli = cells_[i];
        newCenter += (cellCenter[celli] * cellVolume[celli]);
        volume += cellVolume[celli];
    }

    newCenter /= volume;

    meshGeometry_.setCenter(newCenter);

    Info << "Volume avg rotor Center: "<<newCenter<<endl;
}

void rotorMesh::findRotorNormal(rotorGeometry& rotorGeometry)
{

    const vectorField& cellCenter = mesh_.C();

    vector vecAbove = rotorGeometry.direction();


    vector newNormal(Zero);
    
    forAll(cells_, i)
    {   
        label celli = cells_[i];
        vector vecA = cellCenter[celli]-meshGeometry_.center();
    
        forAll(cells_, j)
        {
            label cellj = cells_[j];
            vector vecB = cellCenter[cellj]-meshGeometry_.center();

            vector probe = (vecA ^ vecB);
            if((probe & vecAbove) < 0)
            {
                probe *=-1;
            }
            newNormal += probe;
        }
        
    }

    newNormal.normalise();

    meshGeometry_.setDirection(newNormal);
    rotorGeometry.setDirection(newNormal);

    Info << "Volume avg rotor Direction: "<<newNormal<<endl;
}

void rotorMesh::findRotorRadius()
{
    
    const vectorField& cellCenter = mesh_.C();
    scalar maxRadSqr=0;
    forAll(cells_, i)
    {   
        label celli = cells_[i];
        scalar radsqr = magSqr(cellCenterLC_[i]);
        if(radsqr>maxRadSqr)
        {
            maxRadSqr=radsqr;
        }
    }

    meshGeometry_.setRadius(sqrt(maxRadSqr));

    Info << "Max rotor radius: "<<meshGeometry_.radius()<<endl;
}


}