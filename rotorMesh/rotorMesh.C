#include "rotorMesh.H"

#include "rotorMesh.H"
#include "cellSet.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "unitConversion.H"

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
    rotor_(),
    cellsName_(),
    findClosestCenter_(false),
    correctGeometry_(false)
{

}

void rotorMesh::build(rotorGeometry& rotorGeometry)
{
    this->clear();

    switch (selMode_)
    {
    case selectionMode::smGeometry :

        rotor_= rotorGeometry; //Use as input data

        this->tryUpdateCenter(); 

        this->createMeshSelection(); //Create cell selection from a disk

        this->computeCellsArea();

        //Find "new" geometry --
        this->findRotorCenter();
        this->findRotorRadius();
        this->findRotorNormal();
        this->tryCorrectGeometry(rotorGeometry); //Update data if required

        break;
    case selectionMode::smCellZone :
    case selectionMode::smCellSet :

        this->loadMeshSelection(); //Load mesh from cellSet

        //Find geometry---------
        this->findRotorCenter();
        this->findRotorNormal();
        this->findRotorRadius();

        this->computeCellsArea(); //

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
        correctGeometry_ = true;

        break;
    case selectionMode::smCellSet :

        ok &= dict.readEntry("cellSet",cellsName_);
        Info<< "Reading cellSet from: " << cellsName_<<endl;

        //If rotorMesh is from cellset geometry may not be provided
        correctGeometry_ = true;

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
        vector oldCenter = rotor_.center;
        label centercell = mesh_.findCell(oldCenter);
        if(centercell == -1)
        {
            //Maybe check before if point is inside of boundaries
            Info << "Rotor center is outside of mesh" << endl;
            return;
        }
        rotor_.center = mesh_.C()[centercell];


        Info << "Using rotor center: " << rotor_.center
            << " instead of: " << oldCenter 
            << endl;
    }
}
void rotorMesh::tryCorrectGeometry(rotorGeometry& rotorGeometry)
{
    if(correctGeometry_)
    {
        rotorGeometry.radius = rotor_.radius;
        rotorGeometry.center = rotor_.center;
        rotorGeometry.direction = rotor_.direction;
    }
}
void rotorMesh::createMeshSelection()
{

    const vector& rotorDir = rotor_.direction;
    const vector& rotorCenter = rotor_.center;
    const scalar& radius = rotor_.radius;
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
        vector cel2cent = cellCentroid - rotorCenter;
        scalar distanceSqr = magSqr(cel2cent);

        if( distanceSqr <= radiusSqr )
        {
            cells_.append(celli);
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
    const vector& rotorDir = rotor_.direction;
    const scalar& radius = rotor_.radius;

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

    rotor_.center = newCenter;

    Info << "Volume avg rotor Center: "<<newCenter<<endl;
}

void rotorMesh::findRotorNormal()
{

    const vectorField& cellCenter = mesh_.C();
    const vector vecAbove(1000,1000,1000);

    vector newNormal(Zero);
    
    forAll(cells_, i)
    {   
        label celli = cells_[i];
        vector vecA = cellCenter[celli]-rotor_.center;
    
        forAll(cells_, j)
        {
            label cellj = cells_[j];
            vector vecB = cellCenter[cellj]-rotor_.center;

            vector probe = (vecA ^ vecB);
            if((probe & vecAbove) < 0)
            {
                probe *=-1;
            }
            newNormal += probe;
        }
        
    }

    newNormal.normalise();

    rotor_.direction = newNormal;

    Info << "Volume avg rotor Direction: "<<newNormal<<endl;
}

void rotorMesh::findRotorRadius()
{
    
    const vectorField& cellCenter = mesh_.C();
    scalar maxRadSqr=0;
    forAll(cells_, i)
    {   
        label celli = cells_[i];
        scalar radsqr = magSqr(cellCenter[celli]-rotor_.center);
        if(radsqr>maxRadSqr)
        {
            maxRadSqr=radsqr;
        }
    }


    rotor_.radius = sqrt(maxRadSqr);

    Info << "Max rotor radius: "<<rotor_.radius<<endl;
}


}