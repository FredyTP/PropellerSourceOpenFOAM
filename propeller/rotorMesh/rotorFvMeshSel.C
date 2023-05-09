#include "rotorFvMeshSel.H"
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

defineTypeNameAndDebug(rotorFvMeshSel,0);

const Enum
<
    rotorFvMeshSel::selectionMode
>
rotorFvMeshSel::selectionModeNames_
({
    {selectionMode::smNone, "none"},
    {selectionMode::smGeometry, "geometry"},
    {selectionMode::smCellSet, "cellSet"},
    {selectionMode::smCellZone, "cellZone"},
    {selectionMode::smCylinder, "cylinder"}
});




rotorFvMeshSel::rotorFvMeshSel
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    selMode_(selectionMode::smNone),
    cells_(),
    meshGeometry_(),
    cellsName_(),
    findClosestCenter_(false),
    correctGeometry_(false)
{

}

void rotorFvMeshSel::build(rotorGeometry& rotorGeometry)
{
    Info<<endl;
    Info<<"Building Mesh Selection: "<<endl;
    Info.stream().incrIndent();
    // Fix rotor center specification
    //Provide geometry data if not data specified or ask for correction
    this->clear();

    switch (selMode_)
    {
    case selectionMode::smGeometry :
    case selectionMode::smCylinder :
    {
        meshGeometry_= rotorGeometry; //Use as input data

        const bool rgReady = rotorGeometry.innerRadius().isReady()
                && rotorGeometry.radius().isReady() 
                && rotorGeometry.center().isReady()  
                && rotorGeometry.direction().isReady();

        //Check if geometry is enough to build mesh
        if(!rgReady)
        {
            FatalErrorInFunction
                << "Trying to select rotor mesh with incomplete geometry: "<<endl
                << rotorGeometry
                <<exit(FatalError);
        }

        //Update to closest center if it's asked to (revisar esto)
        this->tryUpdateCenter(); 

        if(selMode_==selectionMode::smGeometry)
        {
            this->createMeshSelection(); //Create cell selection from a disk
        }
        else
        {
            this->createMeshSelectionCilinder();
        }

        //Find "new" geometry --
        this->findRotorRadius();

        if(meshGeometry_.radius().get() < rotorGeometry.radius().get())
        {
            FatalErrorInFunction
                << "Mesh Selection radius is smaller than rotor radius: "
                << "("<<meshGeometry_.radius().get() <<" < "<<rotorGeometry.radius().get()<<")"<<endl
                << exit(FatalError);
        }
        //Dont find rotor center because already selected cells from spec or aprox center
        //this->findRotorCenter();
        //this->findRotorNormal(rotorGeometry);
        this->tryCorrectGeometry(rotorGeometry); //Update data if required


    }
        break;
    case selectionMode::smCellZone :
    case selectionMode::smCellSet :

        //Check for input information
        if(!rotorGeometry.direction().isReady())
        {
            FatalErrorInFunction
            <<"No defined disk normal for rotor mesh selection"
            <<exit(FatalError);
        }

        this->loadMeshSelection(); //Load mesh from cellSet

        //Find geometry---------
        this->findRotorCenter();
        this->findRotorNormal(rotorGeometry); //ROTOR NORMAL ALWAYS UPDATES ROTOR GEOMETRY NORMAL
        this->findRotorRadius();

        //Update to closest center if it's asked to (revisar esto)
        this->tryUpdateCenter();
 
        //For mesh selection, the used radius is the maximum
        //(It is supposed that when a selection is used, a smooth circular disk is provided)
        meshGeometry_.radius().set(maxRadius());
        this->tryCorrectGeometry(rotorGeometry); //Use geometry as output
    break;

    default:
    break;
    }
    
    built_ = true;
    Info.stream().decrIndent();
    

}
void rotorFvMeshSel::clear()
{
    cells_.clear();
    built_ = false;
}
bool rotorFvMeshSel::read(const dictionary &dict)
{
    Info<<endl;
    Info<<"Reading Rotor mesh selection config: "<<endl;
    Info.stream().incrIndent();
    bool ok = true;

    selMode_ = selectionModeNames_.getOrDefault("selectionMode",dict,selectionMode::smGeometry);

    switch (selMode_)
    {
    case selectionMode::smGeometry :
    case selectionMode::smCylinder :
        
        findClosestCenter_ = dict.getOrDefault<bool>("closestCenter",false);
        //includeIfVertex_ = dict.getOrDefault<bool>("includeIfVertex",true);
        indent(Info)<<"- Selection mode: geometry"<<endl;
        indent(Info)<<"- FindClosestCenter: "<< findClosestCenter_<<endl;

        //If rotorFvMeshSel is build from geometry, the resulting mesh geometry
        // may not be used to update the provided geometry
        correctGeometry_ = false;//dict.getOrDefault<bool>("correctGeometry",false);
        indent(Info)<<"- Correct rotor geometry: "<<correctGeometry_<<endl;

        break;

    case selectionMode::smCellZone :
        ok &= dict.readEntry("cellZone",cellsName_);
        findClosestCenter_ = dict.getOrDefault<bool>("closestCenter",false);
        indent(Info)<< "- Selection mode: cellZone(" << cellsName_<<")"<<endl;
        indent(Info)<<"- FindClosestCenter: "<< findClosestCenter_<<endl;

        //If rotorFvMeshSel is from cellset geometry may not be provided
        correctGeometry_ = dict.getOrDefault<bool>("correctGeometry",false);

        indent(Info)<<"- Correct rotor geometry: "<<correctGeometry_<<endl;

        break;
    case selectionMode::smCellSet :

        ok &= dict.readEntry("cellSet",cellsName_);
        findClosestCenter_ = dict.getOrDefault<bool>("closestCenter",false);
        indent(Info)<< "- Selection mode: cellSet(" << cellsName_<<")"<<endl;
        indent(Info)<<"- FindClosestCenter: "<< findClosestCenter_<<endl;

        //If rotorFvMeshSel is from cellset geometry may not be provided
        correctGeometry_ = dict.getOrDefault<bool>("correctGeometry",false);
        indent(Info)<<"- Correct rotor geometry: "<<correctGeometry_<<endl;

        break;
    default:
        break;
    }
    Info.stream().decrIndent();

    return ok;
}
void rotorFvMeshSel::tryUpdateCenter()
{

    if(findClosestCenter_)
    {
        //Use closest cell centroid to center the plane used
        //to find rotor geometry, usually provides better results
        //in an oriented mesh
        vector oldCenter = meshGeometry_.center().get();
        vector newCenter(0,0,0);
        label centercell = mesh_.findCell(oldCenter);
        label found = 0;
        if(centercell != -1)
        {
            newCenter = mesh_.C()[centercell];
            found=1;
        }
        else
        {
            if(!Pstream::parRun())
            {
                Warning<<"Closest center is outside mesh"<<endl;
                return;
            }
        }
        reduce(found,sumOp<label>());
        if(found==0)
        {
            Warning<<"Center cell is outside boundaries, using provided center to cut the cells"<<endl;
            return;
        }
        else if(found>1)
        {
            Warning<<"Multiple cells found, using provided center to cut the cells"<<endl;
            return;
        }
        else
        {
            reduce(newCenter,sumOp<vector>());
            vector n = meshGeometry_.direction().get();
            scalar d = -n.inner(oldCenter);
            scalar dist = newCenter.inner(n)+d;

            vector v0 = newCenter-oldCenter;
            scalar dir = v0.inner(n);
            if(dir<0) dist*=-1;

            newCenter = oldCenter + n*dist;
            meshGeometry_.center().set(newCenter);
        }
        
        indent(Info) << "- Using rotor center: " << meshGeometry_.center().get()
        << " instead of: " << oldCenter 
        << endl;
        
    }
}
void rotorFvMeshSel::tryCorrectGeometry(rotorGeometry& rotorGeometry)
{
    if(findClosestCenter_)
    {
        rotorGeometry.center().set(meshGeometry_.center().get());
    }
    if(correctGeometry_)
    {
        rotorGeometry.radius().set(meshGeometry_.radius().get());
        rotorGeometry.center().set(meshGeometry_.center().get());
        rotorGeometry.direction().set(meshGeometry_.direction().get());
        //Check for no refPsi direction (?)
    }
    else
    {
        //Add missing geometry information
        //Useful when selecting from mesh
        //and some parameters want to be specified like radius
        //when a mesh selection is bigger than the desired radius
        //And others want to be deduced like center or direction

        if(!rotorGeometry.radius().isSet())
        {
            rotorGeometry.radius().set(meshGeometry_.radius().get());
        }
        if(!rotorGeometry.center().isSet())
        {
            rotorGeometry.center().set(meshGeometry_.center().get());
        }
        if(!rotorGeometry.direction().isSet())
        {
            rotorGeometry.direction().set(meshGeometry_.direction().get());
        }
    }
}

void rotorFvMeshSel::syncCellData()
{
    //Join all core cell number
    label allCoreCells = cells_.size();
    reduce(allCoreCells,sumOp<label>());

    //Obtain number of cells for each core
    parNcells.resize(Pstream::nProcs(),0);
    parNcells[Pstream::myProcNo()]=cells_.size();

    reduce(parNcells,sumOp<labelList>());
     //Out total number of cells
    indent(Info)<<"- Total selected cells: "<<allCoreCells<<". In each core: "<<parNcells<<endl;

    //Get cell list for each core
    List<List<label>> coreCells;
    parCells.resize(Pstream::nProcs());
    forAll(parCells,i)
    {
        parCells[i].resize(parNcells[i],0);
    }
    parCells[Pstream::myProcNo()]=cells_;
    forAll(parCells,i)
    {
        reduce(parCells[i],sumOp<labelList>());
    }
}

/**
 * Uncomplete function, need to ensure only 1 mesh layer
*/
void rotorFvMeshSel::createMeshSelectionCilinder()
{

    Warning<<"Cylinder Selection is not currently working"<<endl;

    const vector& rotorDir = meshGeometry_.direction().get();
    const vector& rotorCenter = meshGeometry_.center().get();
    const scalar& radius = meshGeometry_.radius().get();

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
void rotorFvMeshSel::createMeshSelection()
{

    const vector& rotorDir = meshGeometry_.direction().get();
    const vector& rotorCenter = meshGeometry_.center().get();
    const scalar& radius = meshGeometry_.radius().get();
    //TODO: Find better or improve algorithm to select rotor cells
    plane rotorPlane(rotorCenter,rotorDir);

    cuttingPlane cutPlane(rotorPlane,mesh_,false);

    labelList& planeCells = cutPlane.meshCells();

    //User squared radius to compute faster
    scalar radiusSqr = radius * radius;
    //- Local cartesian position
    //Build local coordinate system
    coordSystem::cartesian localCartesianCS(
            meshGeometry_.center().get(), //centerd to local
            meshGeometry_.direction().get(), //z-axis 
            meshGeometry_.psiRef().get()  //x-axis
        );
     
    const labelListList& cellVertex = mesh_.cellPoints(); 
    const pointField& meshPoints = mesh_.points();
    //Select only cells which centroid proyected over the disk results inside the disk
    forAll(planeCells,i)
    { 
        label celli = planeCells[i];

        if(!includeIfVertex_)
        {
            //Include cell if the centroid proyects inside the rotor
            vector cellCentroid = mesh_.C()[celli];
            vector localPos = localCartesianCS.localPosition(cellCentroid);
            localPos.z()=0; //set on rotor plane 
            scalar distanceSqr = magSqr(localPos);

            if( distanceSqr <= radiusSqr )
            {
                cells_.append(celli);
            }
        }
        else
        {
            //Include if any cell vertex lies inside the rotor
            forAll(cellVertex[celli],j)
            {
                vector point = meshPoints[cellVertex[celli][j]];
                vector localPos = localCartesianCS.localPosition(point);
                localPos.z()=0; //set on rotor plane 
                scalar distanceSqr = magSqr(localPos);

                if( distanceSqr <= radiusSqr )
                {
                    cells_.append(celli);
                    break;
                }
            }

        }

    }

    //Sync selection data arround processes
    this->syncCellData();

}
void rotorFvMeshSel::loadMeshSelection()
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

    this->syncCellData();

    
}
/*void rotorFvMeshSel::computeCellsArea()
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

}*/

/*void rotorFvMeshSel::computeCellsAreaVoronoid()
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

}*/

void rotorFvMeshSel::findRotorCenter()
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

    
    reduce(newCenter,sumOp<vector>());
    reduce(volume,sumOp<scalar>());
    newCenter /= volume;
    meshGeometry_.center().set(newCenter);

    indent(Info) << "- Volume average rotor center: "<<newCenter<<endl;
}

void rotorFvMeshSel::findRotorNormal(rotorGeometry& rotorGeometry)
{

    const vectorField& cellCenter = mesh_.C();

    vector vecAbove = rotorGeometry.direction().get();


    vector newNormal(Zero);
    
    forAll(cells_, i)
    {   
        label celli = cells_[i];
        vector vecA = cellCenter[celli]-meshGeometry_.center().get();
    
        forAll(cells_, j)
        {
            label cellj = cells_[j];
            vector vecB = cellCenter[cellj]-meshGeometry_.center().get();

            vector probe = (vecA ^ vecB);
            if((probe & vecAbove) < 0)
            {
                probe *=-1;
            }
            newNormal += probe;
        }
        
    }

    reduce(newNormal,sumOp<vector>());
    newNormal/=mag(newNormal);

    meshGeometry_.direction().set(newNormal);
    rotorGeometry.direction().set(newNormal);

    indent(Info) << "- Volume average rotor direction: "<<newNormal<<endl;
}

void rotorFvMeshSel::findRotorRadius()
{
    //Build local coordinate system
    coordSystem::cartesian localCartesianCS(
            meshGeometry_.center().get(), //centerd to local
            meshGeometry_.direction().get(), //z-axis 
            meshGeometry_.psiRef().get()  //x-axis
        );

    //Cell center data
    const vectorField& cellCenter = mesh_.C();

    //Cell points index
    const labelListList& cellpoints = mesh_.cellPoints();

    //Mesh points data
    const pointField& meshPoints = mesh_.points();

    scalar maxCenterRadiusSqr = 0;
    scalar maxRadSqr=0;
    forAll(cells_, i)
    {   
        label celli = cells_[i];
        const labelList& cellPoints = cellpoints[celli];
        const point& center = cellCenter[celli];

        //Find cell center radius and check if biggest
        vector localP = localCartesianCS.localPosition(center);
        localP.z()=0;
        scalar cradsqr = magSqr(localP);
        if(cradsqr>maxCenterRadiusSqr)
        {
            maxCenterRadiusSqr=cradsqr;
        }

        forAll(cellPoints,j)
        {
            //Find cell vertex radius and check if biggest
            const point& pij = meshPoints[cellPoints[j]];
            scalar radsqr = magSqr(localCartesianCS.localPosition(pij));
            if(radsqr>maxRadSqr)
            {
                maxRadSqr=radsqr;
            }
        }
    }

    maxPointRadius = sqrt(maxRadSqr);

    reduce(maxPointRadius,maxOp<scalar>());
    reduce(maxCenterRadiusSqr,maxOp<scalar>());

    meshGeometry_.radius().set(sqrt(maxCenterRadiusSqr));

    indent(Info) << "- Max vertex radius: "<<maxPointRadius<<endl;
    indent(Info) << "- Max cell center radius: "<<meshGeometry_.radius().get()<<endl;
}


}