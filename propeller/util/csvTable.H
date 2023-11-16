#ifndef PROPELLER_UTIL_CSV_TABLE
#define PROPELLER_UTIL_CSV_TABLE

#include <List.H>
#include <autoPtr.H>
#include <string>


namespace Foam
{

/**
 * @brief Csv table clas for handling communications with csv files
 * 
 * @tparam dataType 
 * @tparam headerType 
 */
template<class dataType, class headerType>
class csvTable
{
private:
    List<List<dataType>> columns;
    List<headerType> header;
    label nColumn=0;

    bool hasHeader = false;
    label skipLines = 0;

    void processLine(std::string& line);
    void processHeader(std::string& line);
    

    dataType convertInput(std::string& val);
    headerType convertHeader(std::string& val);
public:


    /*INSERT*/

    void addRow(List<dataType>& row);
    void setHeader(List<headerType>& newHeader);
    void addCol(List<dataType>& col, headerType headerName = "");

    label nCol() const;
    label nRow() const;
    label index(headerType headerName) const;
    
    bool readFile(fileName path);   

    List<dataType> col(label coli) const ;
    List<dataType> col(headerType headerName) const;

    List<dataType> row(label rowi) const;

    /**
     * Returns a column as list of list for compatibility
    */
    List<List<dataType>> col2(label coli) const;
    List<List<dataType>> col2(headerType headerName) const;

    csvTable(bool hasHeader_=false,label skipLines_ = 0);
    csvTable(const dictionary& dict);
    ~csvTable();

    template<class A, class B>
    friend Ostream& operator<<(Ostream& os, const csvTable<A,B>& table);
    
};

}

#include "csvTable.C"
#endif //PROPELLER_UTIL_CSV_TABLE

