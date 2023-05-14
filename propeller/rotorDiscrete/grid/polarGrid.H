#ifndef PROPELLER_ROTOR_DISCRETE_GRID_POLAR_GRID
#define PROPELLER_ROTOR_DISCRETE_GRID_POLAR_GRID

#include "scalar.H"
#include "List.H"
#include "gridCell.H"
#include "rotorGrid.H"
#include "ijkAddressing.H"
#include "cylindricalCS.H"
namespace Foam
{

class polarGrid : public rotorGrid, public ijkAddressing
{
public:
    
    TypeName("polarGrid");
    
    polarGrid() = default;

    polarGrid
    (
        const dictionary& dict,
        const rotorGeometry& geometry,
        const rotorFvMeshSel &rotorFvMeshSel,
        const bladeModelS& bladeModel,
        scalar nBlades
    );

    //polarGrid(List<scalar>& radius, List<scalar> theta);

    void setRotation(scalar theta0){}
    void assignFvCells();
    void build();

    gridCell& cell(label i, label j) { return cells_[index(i,j,0)];}
    const List<scalar>& radius() const {return radius_;}
    const List<scalar>& theta() const {return theta_;}

private:
    
    List<scalar> radius_;
    List<scalar> theta_;

    //Acces with ijkAdressing
    void buildGrid();

};

}

#endif /* PROPELLER_ROTOR_DISCRETE_GRID_POLAR_GRID */
