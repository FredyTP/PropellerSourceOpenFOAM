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
    if(row.size()!=nColumn && nColumn>0)
    {
        //Invalid row
        return;
    }
    table->push_back(row);
    if(nColumn==0)
    {
        nColumn = row.size();
    }
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
        header->push_back(cell);
    }
    if(nColumn==0)
    {
        nColumn = header->size();
    }
    
}

template <class dataType, class headerType>
label csvTable<dataType, headerType>::nCol()
{
    if(!table.good())
    {
        return 0;
    }
    if(this->nRow()>0)
    {
        return (*table)[0].size();
    }
    else
    {
        return 0;
    }
}

template <class dataType, class headerType>
label csvTable<dataType, headerType>::nRow()
{
    return table.good()? table->size() : 0;
}

template <class dataType, class headerType>
label csvTable<dataType, headerType>::index(headerType headerName)
{
    if(!header.good())
    {
        return -1;
    }
    for(label i = 0 ; i < header->size();i++)
    {
        if((*header)[i]==headerName)
        {
            return i;
        }
    }
    return header->size();
}

template <class dataType, class headerType>
bool csvTable<dataType, headerType>::readFile(fileName path)
{
    using namespace std;

    table = autoPtr<List<List<dataType>>>::New();
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
List<dataType> csvTable<dataType, headerType>::col(label coli)
{
    List<dataType> column(0);
    if(coli>=(*table)[0].size())
    {
        return column;
    }
    column.resize(table->size());
    for(label i = 0 ; i < table->size(); i++)
    {
        column[i] = (*table)[i][coli];
    }

    return column;
}

template <class dataType, class headerType>
List<dataType> csvTable<dataType, headerType>::col(headerType headerName)
{
    return this->col(index(headerName));
}

template <class dataType, class headerType>
List<List<dataType>> csvTable<dataType, headerType>::col2(label coli)
{
    List<List<dataType>> column(0);
    if(coli>(*table)[0].size())
    {
        return column;
    }
    column.resize(table->size());
    for(label i = 0 ; i < table->size(); i++)
    {
        column[i].resize(1);
        column[i][0] = (*table)[i][coli];
    }

    return column;
}

template <class dataType, class headerType>
List<List<dataType>> csvTable<dataType, headerType>::col2(headerType headerName)
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
        std::runtime_error e(val);
        throw e;
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
        std::runtime_error e(val);
        throw e;
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
    if(table.header.good())
    {
        for(label i = 0; i < table.header->size(); i++)
        {
            if(i!=0)
            {
                os<<", ";
            }
            os<<(*table.header)[i];
        }
        os<<endl;
    }
    if(table.table.good())
    {
        for(label i = 0; i< table.table->size();i++)
        {
            const List<A>& row = (*table.table)[i];
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

