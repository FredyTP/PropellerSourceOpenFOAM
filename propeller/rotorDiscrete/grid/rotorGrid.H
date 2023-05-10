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

namespace Foam
{


class rotorGrid
{
public:

    rotorGrid(const coordSystem::cylindrical &cylCS, scalar innerRadius, scalar radius);

    virtual void assignFvCells(const vectorField& cellCenter, const scalarField& weights, const labelList& cellis) = 0;
    virtual void build() = 0;
    void setCenterFromClosestCell(const vectorField& cellCenter);

    const PtrList<gridCell>& cells() const { return cells_;}
    PtrList<gridCell>& cells() { return cells_;}
    const List<vector>& centers() const {return centers_;}

    static tensor bladeLocalFromPoint(const coordSystem::cylindrical &cylCS, const point &localPoint);

    static autoPtr<rotorGrid> New(const dictionary &dict, const coordSystem::cartesian& carCS, const coordSystem::cylindrical &cylCS, scalar innerRadius, scalar radius);

    virtual ~rotorGrid() = default;
protected:
    coordSystem::cylindrical cylCS_;

    scalar minRadius_;
    scalar maxRadius_;


    PtrList<gridCell> cells_;
    List<vector> centers_;
};

}

#endif /* PROPELLER_ROTOR_DISCRETE_ROTOR_GRID */


#endif /* PROPELLER_ROTOR_DISCRETE_GRID_ROTOR_GRID_20COPY */
