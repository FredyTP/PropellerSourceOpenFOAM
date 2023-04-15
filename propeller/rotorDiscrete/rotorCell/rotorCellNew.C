#include "rotorCell.H"
#include "rotorCenterCell.H"
#include "rotorTriCell.H"

namespace Foam
{
const Enum
<
    rotorCell::integrationMode
>
rotorCell::integrationModeNames_
({
    {integrationMode::imCenter, "center"},
    {integrationMode::imTri, "tri"},
    {integrationMode::imTriCenter, "tricenter"}
});

autoPtr<rotorCell> rotorCell::New(integrationMode type, label center, List<label>& vertex, const List<point> &points, label cellIdx)
{
    if(type == integrationMode::imCenter)
    {
        return autoPtr<rotorCell>::NewFrom<rotorCenterCell>(center,vertex,points,cellIdx);
    }
    else if(type == integrationMode::imTri)
    {
        //Triangulated mesh cells without center 
        return autoPtr<rotorCell>::NewFrom<rotorTriCell>(center,vertex,points,cellIdx,false);
    }
    else if(type == integrationMode::imTriCenter)
    {
        //Triangulated mesh cells with center included
        return autoPtr<rotorCell>::NewFrom<rotorTriCell>(center,vertex,points,cellIdx,true);
    }
    else
    {
        FatalErrorInFunction << "Rotor cell type: "
        <<type <<" does not exit, try with: center or tri"
        <<exit(FatalError);
        return autoPtr<rotorCell>();
    }
}

}

