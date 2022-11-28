#include "rotorDiscrete.H"
#include "cellSet.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "unitConversion.H"

namespace Foam
{

defineTypeNameAndDebug(rotorDiscrete,0);

const Enum
<
    rotorDiscrete::selectionMode
>
rotorDiscrete::selectionModeNames_
({
    {selectionMode::smGeometry, "geometry"},
    {selectionMode::smCellSet, "cellSet"},
});




rotorDiscrete::rotorDiscrete(
    scalar radius,
    const fvMesh &mesh,
    const dictionary &dict) 
:
    mesh_(mesh),
    radius_(radius)
{
    this->read(dict);
}

bool rotorDiscrete::read(const dictionary &dict)
{
    Info<<"Reading Rotor geometry"<<endl;

    bool ok = true;

    ok &= selectionModeNames_.readEntry("selectionMode",dict,selMode_);

    switch (selMode_)
    {
    case selectionMode::smGeometry:
        ok &= dict.readEntry("center", rotorCenter_);
        ok &= dict.readEntry("direction", rotorDir_);
        rotorDir_.normalise();
        ok &= dict.readEntry("closestCenter",isClosestCenter_);
          
        if(ok)
        {
            Info<<"Radius: "<< radius_<<endl;
            Info<<"Center: "<< rotorCenter_<<endl;
            Info<<"Direction: "<< rotorDir_<<endl;
            Info<<"FindClosestCenter: "<< isClosestCenter_<<endl;
            bool selectMesh = dict.getOrDefault<bool>("createSelection",true);
            this->updateCenter();
           
            if(selectMesh)
            {
                this->createMeshSelection();
            }
            else{
                Info<<"Not creating mesh selection from geometry, use this for debug purpose"<<endl;
            }
            psiOrigin_ = vector(rotorDir_.y(),-rotorDir_.x(),0);
            psiOrigin_.normalise();
            coordSys_ = coordSystem::cylindrical(realCenter_,rotorDir_,psiOrigin_);
            

            cylCellCenter.resize(cells_.size());

            forAll(cells_,i)
            {
                cylCellCenter[i]=coordSys_.localPosition(mesh_.C()[cells_[i]]);
            }
            
        }
        break;

    case selectionMode::smCellSet:
    //THIS MODE NOT AVAILABLE YET!!
        /* code */
        word zoneName("");

        ok &= dict.readEntry("cellSet",zoneName);
        if(ok)
        {
            cells_ = cellSet(mesh_, zoneName).sortedToc();
        }
        break;
    }

    this->findArea();

    return ok;
}
void rotorDiscrete::updateCenter()
{

    if(isClosestCenter_)
    {
        //Use closest cell centroid to center the plane used
        //to find rotor geometry, usually provides better results
        //in an oriented mesh

        label centercell = mesh_.findCell(rotorCenter_);
        realCenter_ = mesh_.C()[centercell];


        Info << "Using rotor center: " << realCenter_
            << " instead of: " << rotorCenter_ 
            << endl;
    }
    else
    {
        realCenter_ = rotorCenter_;
    }
}
void rotorDiscrete::createMeshSelection()
{

    //TODO: Find better or improve algorithm to select rotor cells
    plane rotorPlane(realCenter_,rotorDir_);

    cuttingPlane cutPlane(rotorPlane,mesh_,false);

    labelList& planeCells = cutPlane.meshCells();

    //User squared radius to compute faster
    scalar radiusSqr = radius_ * radius_;

    //Select only cells which centroid to the rotor center <= radius
    forAll(planeCells,i)
    {
        
        label celli = planeCells[i];

        vector cellCentroid = mesh_.C()[celli];
        vector cel2cent = cellCentroid - realCenter_;
        scalar distanceSqr = magSqr(cel2cent);

        if( distanceSqr <= radiusSqr )
        {
            cells_.append(celli);
        }
    }
    Info<<"Selected cells: "<<cells_.size()<<endl;

}

void rotorDiscrete::findArea()
{
 // NO IDEA !!! :(
    scalar idealArea = constant::mathematical::twoPi * radius_*radius_;

    //set area size
    area_.resize(cells_.size());

    area_= 0.0;
    diskArea_ = 0.0;
    /*const labelListList& cellFaces = mesh_.cellFaces();
    area_.resize(cells_.size());
    diskArea_=0;
    forAll(cells_,i)
    {
        area_[i]=0;
        label celli = cells_[i];
        const labelList faces = cellFaces[celli];

        forAll(faces,j)
        {
            label facej = faces[j];
            vector projArea= mesh_.Sf()[facej] & rotorDir_;
            area_[i] += 0.5*mag(projArea);
        }
        diskArea_+=area_[i];
    }*/

    const label nInternalFaces = mesh_.nInternalFaces();
    const vectorField& Sf = mesh_.Sf();
    const scalarField& magSf = mesh_.magSf();


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
            area_[own] += std::abs(Sf[facei] & rotorDir_);
        }
        else if ((own == -1) && (nbr != -1))
        {
            area_[nbr] += std::abs(Sf[facei] & rotorDir_);
        }
    }

    forAll(area_,i)
    {
        diskArea_+=area_[i];
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
    Info<< "Ideal disk area" << idealArea<<endl;
    Info<< "Disk Area: " << diskArea_ << endl;    

}

}