#include "airfoilModelList.H"
#include "interpolatedAirfoil.H"


Foam::airfoilModelList::airfoilModelList(const dictionary& dict)
: PtrList<airfoilModel>()
{
    wordList airfoilNames(dict.toc());
    Info<<endl;
    Info<< "Creating airfoils:" <<endl;

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

        forAll(airfoilNames,i)
        {
            interpolatedAirfoil* iAirfoil = dynamic_cast<interpolatedAirfoil*>(this->get(i));
            if(iAirfoil!=nullptr)
            {
                iAirfoil->build(*this);
            }
            this->get(i)->exportAirfoil();
        }
    }

}

const Foam::airfoilModel* Foam::airfoilModelList::getAirfoil(const word name) const
{
    forAll(*this,i)
    {
        if(name.compare(this->get(i)->airfoilName())==0)
        {
            return this->get(i);
        }
    }
    return nullptr;
}
const Foam::airfoilModel* Foam::airfoilModelList::getAirfoil(label index) const
{
    return this->get(index);
}

