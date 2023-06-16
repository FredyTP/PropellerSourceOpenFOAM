#ifndef PROPELLER_ROTOR_GRID_MESH_GRID
#define PROPELLER_ROTOR_GRID_MESH_GRID


namespace Foam
{

class meshGrid : public rotorGrid
{
public:
    
    TypeName("meshGrid");
    
    polarGrid() = default;

    meshGrid
    (
        const dictionary& dict,
        const rotorGeometry& geometry,
        const rotorFvMeshSel &rotorFvMeshSel,
        const bladeModelS* bladeModel,
        scalar nBlades
    );

    void assignFvCells();

    const List<scalar>& radius() const {return radius_;}
    const List<scalar>& theta() const {return theta_;}


private:

    scalar area_;
    void buildGrid();

};

}

#endif /* PROPELLER_ROTOR_GRID_MESH_GRID */
