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
#include "dictionary.H"
#include "autoPtr.H"

namespace Foam
{


class rotorGrid
{
public:
    enum sampleMode
    {
        spCenter,
        spClosestCell,
        spCellMean
    };
    static const Enum<sampleMode> sampleModeNames_;

    TypeName("rotorGrid");

    declareRunTimeSelectionTable
    (
        autoPtr,
        rotorGrid,
        dictionary,
        (
            const dictionary& dict,
            const rotorGeometry& geometry,
            const rotorFvMeshSel &rotorFvMeshSel,
            const bladeModelS& bladeModel,
            scalar nBlades
        ),
        (dict,geometry,rotorFvMeshSel,bladeModel,nBlades)
    );
    
    rotorGrid(const dictionary& dict, const rotorGeometry& geometry, const rotorFvMeshSel &rotorFvMeshSel, label nBlades);

    virtual void assignFvCells() = 0;
    virtual void build() = 0;
    virtual void setRotation(const List<scalar>& thetas)=0;

    void setCenterFromClosestCell();
    void updateCenters();
    const PtrList<gridCell>& cells() const { return cells_;}
    const rotorGeometry& geometry() const { return rotorGeometry_;}
    PtrList<gridCell>& cells() { return cells_;}
    sampleMode samplingMode() const {return sampleMode_;}
    static tensor bladeLocalFromPoint(const coordSystem::cylindrical &cylCS, const point &localPoint);
    static autoPtr<rotorGrid> New(const dictionary &dict, const rotorGeometry& geometry, const rotorFvMeshSel &rotorFvMeshSel, const bladeModelS& bladeModel, scalar nBlades);
    label nCells() const { return cells_.size();}
    virtual ~rotorGrid() = default;
protected:
    const rotorGeometry& rotorGeometry_;
    const rotorFvMeshSel& meshSel_;
    label nBlades_;
    scalar minRadius_;
    scalar maxRadius_;

    sampleMode sampleMode_;
    PtrList<gridCell> cells_;
};

}

#endif /* PROPELLER_ROTOR_DISCRETE_ROTOR_GRID */


#endif /* PROPELLER_ROTOR_DISCRETE_GRID_ROTOR_GRID_20COPY */
