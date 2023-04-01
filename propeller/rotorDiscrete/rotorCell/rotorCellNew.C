#include "rotorCell.H"
#include "rotorCenterCell.H"
#include "rotorTriCell.H"

namespace Foam
{
autoPtr<rotorCell> rotorCell::New(word type, label center, List<label>& vertex, const List<point> &points, label cellIdx)
{
    if(type == "center")
    {
        return autoPtr<rotorCell>::NewFrom<rotorCenterCell>(center,vertex,points,cellIdx);
    }
    else if(type == "tri")
    {
        //Triangulated mesh cells without center 
        return autoPtr<rotorCell>::NewFrom<rotorTriCell>(center,vertex,points,cellIdx,false);
    }
    else if(type == "tricenter")
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

