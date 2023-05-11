#ifndef PROPELLER_ROTOR_DISCRETE_GRID_BLADE_GRID
#define PROPELLER_ROTOR_DISCRETE_GRID_BLADE_GRID

#include "rotorGrid.H"
#include "ijkAddressing.H"
#include "cartesianCS.H"

namespace Foam
{

class bladeGrid : public rotorGrid, public ijkAddressing
{
public:

    bladeGrid(const rotorGeometry& geometry, scalar chord, label nBlades, label nRadius, label nChord);
    virtual void assignFvCells(const vectorField& cellCenter, const scalarField& weights, const labelList& cellis);
    virtual void build();

    void buildBlades();
    void updateTheta(scalar theta0);
    void updateThetas(const List<scalar>& thetas);

    void rotateBlades();
    label ijkIndex(label iBlade, label iRadius, label iChord);
private:
    scalar chord_;
    label nRadius_;
    label nChord_;
    
    List<scalar> theta_;

};
}

#endif /* PROPELLER_ROTOR_DISCRETE_GRID_BLADE_GRID */
