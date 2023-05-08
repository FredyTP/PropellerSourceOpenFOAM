#ifndef PROPELLER_ROTOR_DISCRETE_ROTOR_GRID
#define PROPELLER_ROTOR_DISCRETE_ROTOR_GRID

#include "scalar.H"
#include "List.H"
#include "gridCell.H"
#include "ijkAddressing.H"
#include "cylindricalCS.H"
namespace Foam
{


class rotorGrid : public ijkAddressing
{
public:

    rotorGrid() = default;

    /**
     * @brief Construct a new rotor Grid object with a structured grid
     * 
     * @param nRadius 
     * @param nTheta 
     * @param minRadius 
     * @param maxRadius 
     * @param cylCS 
     */
    rotorGrid
    (
        label nRadius,
        label nTheta, 
        scalar minRadius, 
        scalar maxRadius,
        const coordSystem::cylindrical &cylCS
    );
    //rotorGrid(List<scalar>& radius, List<scalar> theta);


    void assignFvCells(const vectorField& cellCenter, const scalarField& weights, const labelList& cellis);
    void build();
    void setCenterFromClosestCell(const vectorField& cellCenter);

    const PtrList<gridCell>& cells() const { return cells_;}
    PtrList<gridCell>& cells() { return cells_;}
    gridCell& cell(label i, label j) { return cells_[index(i,j,0)];}
    const List<scalar>& radius() const {return radius_;}
    const List<scalar>& theta() const {return theta_;}
    const List<vector>& centers() const {return centers_;}

    static tensor bladeLocalFromPoint(const coordSystem::cylindrical &cylCS, const point &localPoint);
private:
    coordSystem::cylindrical cylCS_;
    scalar minRadius_;
    scalar maxRadius_;
    
    List<scalar> radius_;
    List<scalar> theta_;

    //Acces with ijkAdressing
    PtrList<gridCell> cells_;
    List<vector> centers_;

    
    void buildGrid();

};

}

#endif /* PROPELLER_ROTOR_DISCRETE_ROTOR_GRID */
