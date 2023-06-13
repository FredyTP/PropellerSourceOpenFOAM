#include "readHelper.H"
#include "scalarList.H"
#include "unitConversion.H"
#include "LinearInterpolation.H"
#include "cubicSplineInterpolation.H"
namespace Foam
{
namespace util
{

autoPtr<RegularInterpolation<scalar,scalar,1>> NewInterpolationFromDict
(
    const dictionary& dict,
    word x_name,
    word y_name,
    bool convertToRad,
    bool enableCSV,
    const csvTable<scalar,word>* csv
)
{
    List<scalar> xList;
    List<scalar> yList;

    word from = dict.get<word>("from");
    if(from == "csv")
    {
        if(!enableCSV)
        {
            FatalErrorInFunction
                <<"Cannot read from csv for variable: "<<y_name<<endl;
        }
        
        //Get name on csv
        word xCol = dict.get<word>(x_name);
        word yCol = dict.get<word>(y_name);

        if(csv)
        {
            xList = csv->col(xCol);
            yList = csv->col(yCol);
        }
        else
        {
            csvTable<scalar,word> table(dict);
            xList = table.col(xCol);
            yList = table.col(yCol);
        }

    }
    else if(from == "list")
    {
        xList = dict.get<scalarList>(x_name);
        yList = dict.get<scalarList>(y_name);
    }
    else if(from == "constant")
    {
        xList.resize(0);
        yList.resize(1);
        yList[0] = dict.get<scalar>(y_name);
    }
    else
    {
        FatalErrorInFunction
            <<"From: "<<from<< " doesn't exit. Available: (csv list constant)"
            <<exit(FatalError);
    }

    if(convertToRad)
    {
        forAll(yList,i)
        {
            yList[i]=degToRad(yList[i]);
        }
    }

    bool isCubic = dict.getOrDefault<bool>("cubicSpline",false);
    
    if(xList.size() == 0 && yList.size()==1)
    {
        return autoPtr<RegularInterpolation<scalar,scalar,1>>
            ::NewFrom<ConstantInterpolation<scalar,scalar,1>>(yList.first());
    }
  
    if(isCubic && xList.size() < 3)
    {
        Warning
        <<"Cannot use cubic spline if list.size() < 3 using linear interpolation on: "
        <<xList<<endl;
    }
    if(isCubic && xList.size() >= 3)
    {
        return autoPtr<RegularInterpolation<scalar,scalar,1>>
        ::NewFrom<cubicSplineInterpolation>(xList,yList);
    }
    else
    {
        return autoPtr<RegularInterpolation<scalar,scalar,1>>
        ::NewFrom<LinearInterpolation<scalar,scalar,1>>(FixedList<scalarList,1>({xList}),yList);
    }

}




}
}



