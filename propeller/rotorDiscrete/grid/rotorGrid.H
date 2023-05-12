#ifndef PROPELLER_ROTOR_DISCRETE_GRID_ROTOR_GRID_20COPY
#define PROPELLER_ROTOR_DISCRETE_GRID_ROTOR_GRID_20COPY
#ifndef PROPELLER_ROTOR_DISCRETE_ROTOR_GRID
#define PROPELLER_ROTOR_DISCRETE_ROTOR_GRID

#include "scalar.H"
#include "List.H"
#include "gridCell.H"
#include "ijkAddressing.H"
#include "cylindricalCS.H"
#include "cartesianCS.H"
#include "rotorGeometry.H"
#include "rotorFvMeshSel.H"
#include "bladeModelS.H"

namespace Foam
{


class rotorGrid
{
public:

    rotorGrid(const rotorGeometry& geometry, const rotorFvMeshSel &rotorFvMeshSel, label nBlades);

    virtual void assignFvCells() = 0;
    virtual void build() = 0;
    void setCenterFromClosestCell();

    const PtrList<gridCell>& cells() const { return cells_;}
    PtrList<gridCell>& cells() { return cells_;}
    const List<vector>& centers() const {return centers_;}
    
    static tensor bladeLocalFromPoint(const coordSystem::cylindrical &cylCS, const point &localPoint);

    static autoPtr<rotorGrid> New(const dictionary &dict, const rotorGeometry& geometry, const rotorFvMeshSel &rotorFvMeshSel, const bladeModelS& bladeModel, scalar nBlades);

    virtual ~rotorGrid() = default;
protected:
    const rotorFvMeshSel& meshSel_;
    const rotorGeometry& rotorGeometry_;
    label nBlades_;
    scalar minRadius_;
    scalar maxRadius_;


    PtrList<gridCell> cells_;
    List<vector> centers_;
};

}

#endif /* PROPELLER_ROTOR_DISCRETE_ROTOR_GRID */


#endif /* PROPELLER_ROTOR_DISCRETE_GRID_ROTOR_GRID_20COPY */
