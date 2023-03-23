#ifndef PROPELLER_UTIL_CSV_TABLE
#define PROPELLER_UTIL_CSV_TABLE

#include <List.H>
#include <autoPtr.H>
#include <string>


namespace Foam
{

/**
 * ADD ERROR CHECKING 
*/
template<class dataType, class headerType>
class csvTable
{
private:
    autoPtr<List<List<dataType>>> table;
    autoPtr<List<headerType>> header;
    label nColumn=0;

    bool hasHeader = false;
    label skipLines = 0;

    void processLine(std::string& line);
    void processHeader(std::string& line);
    

    dataType convertInput(std::string& val);
    headerType convertHeader(std::string& val);
public:

    label nCol();
    label nRow();
    label index(headerType headerName);
    
    bool readFile(fileName path);   

    List<dataType> col(label coli);
    List<dataType> col(headerType headerName);

    /**
     * Returns a column as list of list for compatibility
    */
    List<List<dataType>> col2(label coli);
    List<List<dataType>> col2(headerType headerName);

    csvTable(bool hasHeader_=false,label skipLines_ = 0);
    ~csvTable();

    template<class A, class B>
    friend Ostream& operator<<(Ostream& os, const csvTable<A,B>& table);
    
};
}

#include "csvTable.C"
#endif //PROPELLER_UTIL_CSV_TABLE

