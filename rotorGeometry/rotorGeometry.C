#include "rotorGeometry.H"
#include "cellSet.H"


namespace Foam
{

defineTypeNameAndDebug(rotorGeometry,0);

const Enum
<
    rotorGeometry::selectionMode
>
rotorGeometry::selectionModeNames_
({
    {selectionMode::smGeometry, "geometry"},
    {selectionMode::smCellSet, "cellSet"},
});




rotorGeometry::rotorGeometry(
    const fvMesh &mesh,
    const dictionary &dict) :
    mesh_(mesh)
{
    this->read(dict);
}

bool rotorGeometry::read(const dictionary &dict)
{
    bool ok = true;
    const auto& coeffs = dict.subDict("rotorGeometry");

    ok &= selectionModeNames_.readEntry("selectionMode",coeffs,selMode_);

    switch (selMode_)
    {
    case selectionMode::smGeometry:

        ok &= coeffs.readEntry("radius", radius_);
        ok &= coeffs.readEntry("center", rotorCenter_);
        ok &= coeffs.readEntry("direction", rotorDir_);
        ok &= coeffs.readEntry("closestCenter",isClosestCenter_);
          
        if(ok)
        {
            Info<<"Radius: "<< radius_<<endl;
            Info<<"Center: "<< rotorCenter_<<endl;
            Info<<"Direction: "<< rotorDir_<<endl;
            Info<<"FindClosestCenter: "<< isClosestCenter_<<endl;

            this->updateCenter();
            this->createMeshSelection();
        }
        break;

    case selectionMode::smCellSet:
        /* code */
        word zoneName("");

        ok &= coeffs.readEntry("cellSet",zoneName);
        if(ok)
        {
            cells_ = cellSet(mesh_, zoneName).sortedToc();
        }

        break;
    
    }



    return ok;
}
void rotorGeometry::updateCenter()
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
void rotorGeometry::createMeshSelection()
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

        if(distanceSqr<=radiusSqr)
        {
            cells_.append(celli);
        }
    }
    


}

}