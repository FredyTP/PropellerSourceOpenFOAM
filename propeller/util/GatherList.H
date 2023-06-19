#ifndef PROPELLER_UTIL_GATHER_LIST
#define PROPELLER_UTIL_GATHER_LIST

namespace Foam
{
namespace util
{

/**
 * @brief Gathers in single list all cores list
 * 
 * @tparam type 
 * @param list all core list together
 * @param idx idx where each core begins
 */
template<class type>
List<type> GatherList(const List<type>& list, List<label>& idx)
{
    // Gather ncell information
    List<label> sizes(Pstream::nProcs(), 0);
    idx.resize(Pstream::nProcs());
    sizes[Pstream::myProcNo()] = list.size();
    reduce(sizes, sumOp<labelList>());
    label len = sum(sizes);
    label sum = 0;
    forAll(idx, i)
    {
        idx[i] = sum;
        sum += sizes[i];
    }
    List<type> gatheredList(len,Zero);

    auto it = list.begin();
    auto dstIt = gatheredList.begin() + idx[Pstream::myProcNo()];
    while (it != list.end())
    {
        (*dstIt) = (*it);
        it++;
        dstIt++;
    }
    reduce(gatheredList, sumOp<List<type>>());
    return gatheredList;
}
/**
 * @brief Gathers in single list of lists all cores list of lists
 * 
 * @tparam type 
 * @param list all core list of listtogether
 * @param idx idx where each core begins
 */
template<class type>
List<List<type>> GatherListList(const List<List<type>>& list, List<label>& idx)
{
    // Gather ncell information
    List<label> sizeList(list.size(),0);
    forAll(list,i)
    {
        sizeList[i]=list[i].size();
    }

    List<label> allSize = GatherList<label>(sizeList,idx);
    label len = sum(allSize);
    List<List<type>> gatheredList(allSize.size());
    forAll(gatheredList,i)
    {
        gatheredList[i].resize(allSize[i],Zero);
    }


    auto it = list.begin();
    auto dstIt = gatheredList.begin() + idx[Pstream::myProcNo()];
    while (it != list.end())
    {
        (*dstIt) = (*it);
        it++;
        dstIt++;
    }
    
    forAll(gatheredList,i)
    {
        reduce(gatheredList[i], sumOp<List<type>>());
    }
    return gatheredList;
}

}
}
#endif /* PROPELLER_UTIL_GATHER_LIST */
