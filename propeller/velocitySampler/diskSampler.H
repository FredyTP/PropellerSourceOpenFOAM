#ifndef PROPELLER_VELOCITY_SAMPLER_DISK_SAMPLER
#define PROPELLER_VELOCITY_SAMPLER_DISK_SAMPLER


#include "dictionary.H"
#include "autoPtr.H"
#include "rotorGrid.H"
#include "rotorFvMeshSel.H"
#include "runTimeSelectionTables.H"
namespace Foam
{

/**
 * @brief Velocity sample base class, used to provide velocity data to
 * propeller model
 * 
 * 
 * @example
 * 
 * Example of dicts provided to this class
 * 
 * -Fixed velocity sampler
 * 
 *   diskSampler
 *   {
 *   
 *       type            fixedVelocity;
 *       normal          true;
 *       velocity        10;
 *   }
 */
template<class fType>
class diskSampler
{
public:
    TypeName("diskSampler");
    declareRunTimeSelectionTable
    (
        autoPtr,
        diskSampler,
        dictionary,
        (
            const dictionary& dict,
            const rotorGrid* rGrid,
            const rotorFvMeshSel* rMesh
        ),
        (dict,rGrid,rMesh)
    );
    
    virtual ~diskSampler() = default;

    /**
     * @brief Return sampled velocities on rotor discrete positions(rotor positions goes from 0 to N)
     * Do not get confused by celli which usually means fvmesh cell index
     * 
     * @param U fv velocity field
     * @return const vectorField& 
     */
    virtual const Field<fType>& sampleField(const GeometricField<fType, fvPatchField, volMesh>& U) = 0;

    /**
     * @brief Write field of sampled cells used to interpolate velocity field
     * 
     * @param name Propeller's name for output, adding ":sampledCells" to that name.
     *  Ex: propeller1:sampledCells
     */
    virtual void writeSampled(const word& name);
    
    /**
     * @brief Builds a diskSampler from its dictionary, rotorGrid data and fv Mesh selection
     * 
     * @param dict velocity sampler dict
     * @param rDiscrete_ rotor discretization
     * @param rMesh_ fvMesh selection for rotor
     * @return autoPtr<diskSampler>: safe pointer to created sampler
     */
    static autoPtr<diskSampler> New(const dictionary& dict,const rotorGrid* rGrid,const rotorFvMeshSel* rMesh_);

protected:

    /**
     * @brief Construct a new velocity Sampler object, only available to build from derived clases
     * 
     * @param rGrid 
     * @param rMesh
     */
    diskSampler(const rotorGrid* rGrid,const rotorFvMeshSel* rMesh);

    const rotorGrid* rGrid_;
    const rotorFvMeshSel* rMesh_;

    /**
     * @brief velocity sampled vector, to return by reference
     * or be more efficient when constant velocity is required
     * 
     */
    Field<fType> sampledField_;
};


}


#include "diskSampler.H"

#endif /* PROPELLER_VELOCITY_SAMPLER_DISK_SAMPLER */


