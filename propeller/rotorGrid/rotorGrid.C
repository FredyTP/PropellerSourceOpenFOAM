#include "rotorGrid.H"
#include "polarGrid.H"
#include "bladeGrid.H"
#include "mathematicalConstants.H"

#include "runTimeSelectionTables.H"



namespace Foam
{

defineTypeNameAndDebug(rotorGrid,0);
defineRunTimeSelectionTable(rotorGrid, dictionary);

const Enum
<
    rotorGrid::sampleMode
>
rotorGrid::sampleModeNames_
({
    {sampleMode::spCenter, "center"},
    {sampleMode::spClosestCell, "closestCell"},
    {sampleMode::spCellMean,"cellMean"}
});

rotorGrid::rotorGrid(const dictionary &dict, const rotorGeometry &geometry, const rotorFvMeshSel &rotorFvMeshSel, label nBlades)
    :rotorGeometry_(geometry), meshSel_(rotorFvMeshSel)
{
    sampleMode_ = sampleModeNames_.getOrDefault("sampleMode",dict,sampleMode::spClosestCell);
    nBlades_ = nBlades;
    minRadius_ = geometry.innerRadius();
    maxRadius_ = geometry.radius();
}

void rotorGrid::setCenterFromClosestCell()
{
    const auto& cellCenter = meshSel_.mesh().C();
    forAll(cells_, i)
    {
        cells_[i].centerFromClosestCell(cellCenter);
    }
}
void rotorGrid::build()
{
    forAll(cells_,i)
    {
        cells_[i].checkCells(!meshSel_.isSingleCore(), meshSel_.coreNo());
        cells_[i].buildWeigths(!meshSel_.isSingleCore());
    }
    updateCenters();
}
void Foam::rotorGrid::updateCenters()
{
    if(samplingMode() == spCenter || samplingMode() == spCellMean)
    {
        forAll(cells_, i)
        {
            vector center = cells_[i].getCellCenter();
            cells_[i].setCenter(center);
        }
    }
    else if(samplingMode() == spClosestCell)
    {
        setCenterFromClosestCell();
    }

}

tensor rotorGrid::bladeLocalFromPoint(const coordSystem::cylindrical &cylCS, const point &localPoint) 
{
    // z- up, y -outwards from center, x perpendicular y,z (leading edge to trailing edge)
    point global, origin;

    global = cylCS.globalPosition(localPoint);
    origin = cylCS.origin();

    tensor rotTensor(cylCS.R());

    // z-axis
    const vector ax3 = rotTensor.col<2>(); // == e3 (already normalized)

    // y-axis (radial direction)
    vector ax2(global - origin);

    ax2.removeCollinear(ax3);

    const scalar magAxis2(mag(ax2));

    // Trap zero size and colinearity
    if (magAxis2 < SMALL)
    {
        return rotTensor;
    }

    ax2 /= magAxis2; // normalise

    // Replace with updated local axes

    rotTensor.col<0>(ax2 ^ ax3);
    rotTensor.col<1>(ax2);

    return rotTensor;
}

autoPtr<rotorGrid> Foam::rotorGrid::New(const dictionary &dict, const rotorGeometry &geometry, const rotorFvMeshSel &rotorFvMeshSel, const bladeModelS* bladeModel, scalar nBlades)
{
     //Get model Type name (Ex: simpleAirfoil) 
    //From typeNkey from dictionary (airfoilModel)
    const word modelType(dict.get<word>("type")); 

    Info<< "    - Loading " << modelType << endl;

    //Find class contructor in tables
    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            typeName,
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<rotorGrid>(ctorPtr(dict,geometry,rotorFvMeshSel,bladeModel,nBlades));
}

}

