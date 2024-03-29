#ifndef PROPELLER_ROTOR_DISCRETE_GRID_BLADE_GRID
#define PROPELLER_ROTOR_DISCRETE_GRID_BLADE_GRID

#include "rotorGrid.H"
#include "ijkAddressing.H"
#include "cartesianCS.H"
#include "bladeModelS.H"

namespace Foam
{

/**
 * @brief Blade grid for transient simulations
 * 
 */
class bladeGrid : public rotorGrid, public ijkAddressing
{
public:
    TypeName("bladeGrid");
    
    bladeGrid
    (
        const dictionary& dict,
        const rotorGeometry& geometry,
        const rotorFvMeshSel &rotorFvMeshSel,
        const bladeModelS* bladeModel,
        scalar nBlades
    );
    
    virtual void assignFvCells();
    
    void setRotation(const List<scalar>& thetas);

    label ijkIndex(label iBlade, label iRadius, label iChord);

    void writePythonPlotter(word outputName);
    List<scalar> getInitialAzimuth() const;
private:

    void updateTheta(scalar theta0);
    void updateThetas(const List<scalar>& thetas);
    void rotateBlades();
    void buildBladesConstantChord(scalar chord);
    void buildBladesFromBladeModel(const bladeModelS& bladeModel);

    label nRadius_;
    label nChord_;
    
    List<scalar> theta_;

    void checkDistribution();
    std::function<List<scalar>(label)> chordWiseDistribution;

    static std::function<List<scalar>(label)> constantDistribution();
};
}

#endif /* PROPELLER_ROTOR_DISCRETE_GRID_BLADE_GRID */
