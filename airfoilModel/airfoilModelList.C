#include "airfoilModelList.H"



Foam::airfoilModelList::airfoilModelList(const dictionary& dict)
: PtrList<airfoilModel>()
{
    wordList airfoilNames(dict.toc());
    Info<< "Reading airfoils" <<endl;

    if(airfoilNames.size())
    {
        this->setSize(airfoilNames.size());

        forAll(airfoilNames,i)
        {
            const word& airfoilName = airfoilNames[i];
            this->set(
                i,
                airfoilModel::New(airfoilName, dict.subDict(airfoilName))
            );
        }
    }

}

