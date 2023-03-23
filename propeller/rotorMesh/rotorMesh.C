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

        meshGeometry_= rotorGeometry; //Use as input data

        localCartesianCS_ = coordSystem::cartesian
            (
                meshGeometry_.center(), //centerd to local
                meshGeometry_.direction(), //z-axis 
                meshGeometry_.psiRef()  //x-axis
            );
        this->tryUpdateCenter(); 

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

        break;
    case selectionMode::smCellZone :
    case selectionMode::smCellSet :

        this->loadMeshSelection(); //Load mesh from cellSet

        //Find geometry---------
        this->findRotorCenter();
        this->findRotorNormal(rotorGeometry);
        this->findRotorRadius();

        localCartesianCS_ = coordSystem::cartesian
            (
                meshGeometry_.center(), //centerd to local
                meshGeometry_.direction(), //z-axis 
                meshGeometry_.psiRef()  //x-axis
            );

        this->computeLocalPositions();

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
            Info << "Rotor center is outside of mesh" << endl;
            return;
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
        localPos.z()=0; //set on rotor plane
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

    List<label> idx(centers.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(),idx.end(),[&centers](label a, label b)
    {
        return centers[a].x()<centers[b].x();
    });
    
    List<point> sortedCenter(centers.size());
    forAll(centers,i)
    {
        sortedCenter[i]=centers[idx[i]];
    }
    delaunayTriangulation::Voronoid
    (
        sortedCenter,
        vertex,
        voroCells,
        delaunayTriangulation::circularRegion(radius),
        delaunayTriangulation::intersectCircle(radius)
    );

    forAll(voroCells,i)
    {
        auto& poli = voroCells[i];
        if(poli.size()>2)
        {
            scalar area = 0.0;

            label j = poli.size() - 1;
            for (label i = 0; i < poli.size(); i++)
            {
                area += (vertex[poli[j]].x() + vertex[poli[i]].x()) * (vertex[poli[j]].y() - vertex[poli[i]].y());
                j = i;
            }
            area_[idx[i]] = std::abs(area / 2.0);
           
        }
        else
        {
            Info<<"ERROR BUILDING VORONOID MESH"<<endl;
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

    vector vecAbove(1000,1000,1000);

    //If direction is a valid vector
    if(rotorGeometry.direction() != vector())
    {
        vecAbove=rotorGeometry.direction();
    }

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

    Info << "Volume avg rotor Direction: "<<newNormal<<endl;
}

void rotorMesh::findRotorRadius()
{
    
    const vectorField& cellCenter = mesh_.C();
    scalar maxRadSqr=0;
    forAll(cells_, i)
    {   
        label celli = cells_[i];
        scalar radsqr = magSqr(cellCenter[celli]-meshGeometry_.center());
        if(radsqr>maxRadSqr)
        {
            maxRadSqr=radsqr;
        }
    }

    meshGeometry_.setRadius(sqrt(maxRadSqr));

    Info << "Max rotor radius: "<<meshGeometry_.radius()<<endl;
}


}