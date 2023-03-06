#ifndef PROPELLER_SOURCE_H
#define PROPELLER_SOURCE_H

#include "fvOption.H"
#include "fvMesh.H"
#include "dictionary.H"
#include "Time.H"

#include "cylindricalCS.H"
#include "Enum.H"
#include "volFieldsFwd.H"

#include "propellerModel.H"
#include "airfoilModel.H"
#include "airfoilModelList.H"

#include "irregularInterpolation.H"
#include "bladeModelDev.H"
#include "rotorMesh.H"
#include "rotorGeometry.H"

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class propellerSource Declaration
\*---------------------------------------------------------------------------*/

class propellerSource
:
    public fv::option
{
public:

    //Public Enumerations

    //INFO: each enum must include a Enum class to map
    //each enum item to a string, so it can be selected from 
    //IOdicts
protected:

    //Propeller model used to compute source terms
    autoPtr<propellerModel> propellerModel_;
    rotorMesh rotorMesh_;

    //Disk Geometry
    rotorGeometry rotorGeometry_;
public:
    //Defines the runtime name of this class (IOpurposes?)
    TypeName("propeller");

    // Constructors

        //- Construct from components
        propellerSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        propellerSource(const propellerSource&) = delete;

        //- No copy assignment
        void operator=(const propellerSource&) = delete;


    //- Destructor
    virtual ~propellerSource() = default;
    

    //IO
        //- Read source dictionary
        virtual bool read(const dictionary& dict);


    // Source term addition

        //- Add explicit contribution to momentum equation
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to compressible momentum equation
        /*virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );*/

};



}
}

































#endif //PROPELLER_SOURCE