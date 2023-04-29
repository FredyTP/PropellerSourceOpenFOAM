#include "csvTable.H"
#include <fstream>

namespace Foam
{
template <class dataType, class headerType>
void csvTable<dataType, headerType>::processLine(std::string &line)
{
    using namespace std;

    //0 len list
    List<dataType> row(0);

    stringstream ss(line);
    string value;
    while(getline(ss,value,','))
    {
        dataType cell = convertInput(value);
        row.push_back(cell);
    }
    if(row.size() == 0 )
    {
        //Empty line is okey but dont add it
        return;
    }
    if(row.size()!=columns.size() && columns.size()>0)
    {   
        throw (std::runtime_error("invalid row: '"+value+"'"));
    }
 
    this->addRow(row);
}

template <class dataType, class headerType>
void csvTable<dataType, headerType>::processHeader(std::string &line)
{
    header = autoPtr<List<headerType>>::New(0);
    using namespace std;

    stringstream ss(line);
    string value;
    while(getline(ss,value,','))
    {
        headerType cell = convertHeader(value);
        header.append(cell);
    }
    if(columns.size()==0)
    {
        columns.resize(header.size());
    }
    
}

template <class dataType, class headerType>
void csvTable<dataType, headerType>::addRow(List<dataType> &row)
{
    //Check sizes 
    if(columns.size()==0)
    {
        columns.resize(row.size());
    }
    if(columns.size() != row.size())
    {
        throw(std::runtime_error("Invalid row size : "+row.size()));
    }

    for(label i = 0; i < columns.size(); i++)
    {
        columns[i].append(row[i]);
    }
}

template <class dataType, class headerType>
void csvTable<dataType, headerType>::setHeader(List<headerType> &newHeader)
{
    if(hasHeader && (columns.size()==header.size() || columns.size()==0))
    {
        header=newHeader;
    }
}

template <class dataType, class headerType>
void csvTable<dataType, headerType>::addCol(List<dataType> &col, headerType headerName)
{
    if(columns.size()==0)
    {
        if(headerName!="" && hasHeader)
        {
            header.append(headerName);
        }
        columns.append(col);
    }
    else if(col.size() == this->nRow())
    {
        if(headerName!="" && hasHeader)
        {
            header.append(headerName);
        }
        
        columns.append(col);
    }
}

template <class dataType, class headerType>
label csvTable<dataType, headerType>::nRow() const
{
    if(this->nCol()>0)
    {
        return columns[0].size();
    }
    else
    {
        return 0;
    }
}

template <class dataType, class headerType>
label csvTable<dataType, headerType>::nCol() const
{
    return columns.size();
}

template <class dataType, class headerType>
label csvTable<dataType, headerType>::index(headerType headerName) const
{
    return header.find(headerName);
}

template <class dataType, class headerType>
bool csvTable<dataType, headerType>::readFile(fileName path) 
{
    using namespace std;

    columns.clear();
    string line;
    ifstream fileStream(path);
    if(!fileStream.is_open())
    {
        //ERROR READING FILE
        return false;
    }
    
    //Skip lines needed
    if(skipLines>0)
    {
        label linecount = 0;
        while(getline(fileStream,line))
        {
            if(linecount==skipLines)
            {
                break;
            }
        }
    }
    //Read header
    if(hasHeader)
    {
        getline(fileStream,line);
        this->processHeader(line);
    }

    //read document
    while(getline(fileStream,line))
    {
        this->processLine(line);
    }
    return true;
}

template <class dataType, class headerType>
List<dataType> csvTable<dataType, headerType>::col(label coli) const
{
    if(coli < 0 || coli > this->nCol() -1)
    {
        //Index out of bound return empty column
        return List<dataType>();
    }
    return columns[coli];
}

template <class dataType, class headerType>
List<dataType> csvTable<dataType, headerType>::col(headerType headerName) const
{
    return this->col(index(headerName));
}

template <class dataType, class headerType>
List<dataType> csvTable<dataType, headerType>::row(label rowi) const
{
    if(rowi < 0 || rowi > this->nRow() -1)
    {
        //Index out of bound return empty column
        return List<dataType>();
    }
    List<dataType> trow(this->nCol());
    for(label i = 0; i<trow.size();i++)
    {
        trow[i]=columns[i][rowi];
    }
    return trow;
}

template <class dataType, class headerType>
List<List<dataType>> csvTable<dataType, headerType>::col2(label coli) const
{
    List<List<dataType>> column(0);
    if(coli < 0 || coli > this->nCol() -1)
    {
        //Index out of bound return empty column
        return column;
    }

    column.resize(this->nRow());
    for(label i = 0 ; i < this->nRow(); i++)
    {
        column[i].resize(1);
        column[i][0] = (columns)[coli][i];
    }

    return column;
}

template <class dataType, class headerType>
List<List<dataType>> csvTable<dataType, headerType>::col2(headerType headerName) const
{
    return this->col2(index(headerName));
}

template <class dataType, class headerType>
dataType csvTable<dataType, headerType>::convertInput(std::string &val)
{
    dataType valor;

    std::stringstream stream(val);
    stream >> valor;
    if (stream.fail()) {
        throw std::runtime_error("invalid type conversion of value: "+val);
    }
    return valor;
}

template <class dataType, class headerType>
headerType csvTable<dataType, headerType>::convertHeader(std::string &val)
{
    headerType valor;

    std::stringstream stream(val);
    stream >> valor;
    if (stream.fail()) {
        throw std::runtime_error("invalid type conversion of header: "+val);
    }
    return valor;
}

template <class dataType, class headerType>
csvTable<dataType, headerType>::csvTable(bool hasHeader_,label skipLines_)
 : hasHeader(hasHeader_),skipLines(skipLines_)
{

}

template <class dataType, class headerType>
csvTable<dataType, headerType>::~csvTable()
{
    
}
template<class A, class B>
Ostream& operator<<(Ostream& os, const csvTable<A, B> &table)
{
    if(table.header.size()>0)
    {
        for(label i = 0; i < table.header.size(); i++)
        {
            if(i!=0)
            {
                os<<", ";
            }
            os<<table.header[i];
        }
        os<<endl;
    }
    if(table.columns.size()>0)
    {
        for(label i = 0; i< table.nRow();i++)
        {
            const List<A> row = table.row(i);
            for(label j =0; j< row.size();j++)
            {
                if(j!=0)
                {
                    os<<", ";
                }
                os<<(row[j]);
            }
            os<<endl;
        }
    }
    return os;
}

}

