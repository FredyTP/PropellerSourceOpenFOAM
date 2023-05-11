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

    polarGrid() = default;

    /**
     * @brief Construct a new rotor Grid object with a structured grid
     * 
     * @param nRadius 
     * @param nTheta 
     * @param minRadius 
     * @param maxRadius 
     * @param cylCS 
     */
    polarGrid
    (
        const rotorGeometry& geometry,
        label nBlades,
        label nRadius,
        label nTheta
    );
    //polarGrid(List<scalar>& radius, List<scalar> theta);


    void assignFvCells(const vectorField& cellCenter, const scalarField& weights, const labelList& cellis);
    void build();

    gridCell& cell(label i, label j) { return cells_[index(i,j,0)];}
    const List<scalar>& radius() const {return radius_;}
    const List<scalar>& theta() const {return theta_;}

    static tensor bladeLocalFromPoint(const coordSystem::cylindrical &cylCS, const point &localPoint);
private:
    
    List<scalar> radius_;
    List<scalar> theta_;

    //Acces with ijkAdressing
    void buildGrid();

};

}

#endif /* PROPELLER_ROTOR_DISCRETE_GRID_POLAR_GRID */
