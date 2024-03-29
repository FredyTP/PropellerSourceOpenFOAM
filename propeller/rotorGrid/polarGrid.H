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

/**
 * @brief Polar grid type with r-theta cells
 * 
 */
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
        const bladeModelS* bladeModel,
        scalar nBlades
    );

    //polarGrid(List<scalar>& radius, List<scalar> theta);
    void assignFvCells();

    gridCell& cell(label i, label j) { return cells_[index(i,j,0)];}
    const List<scalar>& radius() const {return radius_;}
    const List<scalar>& theta() const {return theta_;}

    label nRadial() const {return radius().size()-1;}
    label nAzimutal() const {return theta().size()-1;}
    void writePythonPlotter(word outputName);
private:
    
    List<scalar> radius_;
    List<scalar> theta_;

    //Acces with ijkAdressing
    void buildGrid();

};

}

#endif /* PROPELLER_ROTOR_DISCRETE_GRID_POLAR_GRID */
