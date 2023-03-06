#ifndef PROPELLER_INTERPOLATION_REGULAR_LINEAR_ITERATORS
#define PROPELLER_INTERPOLATION_REGULAR_LINEAR_ITERATORS
/*class tableIterator
{
public:
    virtual tableIterator next() = 0;
    virtual tableIterator begin() = 0;
    virtual tableIterator end() = 0;

    virtual scalar operator*() = 0; //returns accesed member
    virtual tableIterator& operator++() = 0; //increments position by 1
    virtual bool operator!=(tableIterator& lhs, tableIterator& rhs) = 0; // compares 2 table iterators
};

template<class itType_>
    static label FindIndex(scalar input,itType_ inputs, label& i1, label& i2)
    {

        auto it = inputs.begin();
        auto next = it.next();
        if(input<*it)
        {
            //indexes are both set to closest index
            i1=it.index();
            i2=it.index();
            //out of bounds
            return 0;
        }

        for(;next!=inputs.end();++it)
        {
            next = it.next();
            //If (luckily) there is coincidence
            if(input == *it)
            {
                i1=it.index();
                i2=it.index();
                return 1;
            }
            //Between i and i+1
            if(*it < input && input< *next)
            {
                i1=it.index();
                i2=next.index();
                return 1;
            }
        }

        //Out of bonds
        i1=next.index();
        i2=i1;
        return 0;
    }
*/

class listIterator
{
private:
    const List<scalar>* data;
    label idx;
public:
    listIterator(const List<scalar>* data_) :data(data_)
    {
        idx=0;
    }
    listIterator(const List<scalar>* data_,label beginIdx) :data(data_)
    {
        idx=beginIdx;
    }
    listIterator next()
    {
        return listIterator(data,idx+1);
    }
    listIterator begin()
    {
        return listIterator(data,0);
    }
    listIterator end()
    {
        return listIterator(data,data->size()-1);
    }
    void operator=(const listIterator& other)
    {
        this->data = other.data;
        this->idx = other.idx;
    }
    scalar operator*()
    {
        return (*data)[idx];
    }
    listIterator& operator++()
    {
        idx++;
        return *this;
    }
    scalar index() {return idx;}
    friend bool operator==(listIterator lhs, listIterator rhs) 
    {
        return ((lhs.data) == (rhs.data)) && (lhs.idx == rhs.idx);
    }
    friend bool operator!=(listIterator lhs, listIterator rhs)  
    {
        return !(lhs==rhs);
    }
    
};

class tableIterator2d
{
private:
    const List<List<scalar>>* data;
    label idx1;
    label idx2;
    bool traverse1 = true;
public:
    tableIterator2d(const List<Listscalar>>* data_,label idx1_, label idx2_, bool traverse1_) :data(data_)
    {
        idx1=idx1_;
        idx2=idx2_;
    }
    tableIterator2d(const List<Listscalar>>* data_, bool traverse1_) :data(data_)
    {
        idx1=0;
        idx2=0;
    }

    tableIterator2d next()
    {
        if(traverse1)
        {
            return tableIterator2d(data,idx1+1,idx2,traverse1);
        }
        else{
            
            return tableIterator2d(data,idx1,idx2+1,traverse1);
        }
    } 
    tableIterator2d begin()
    {
        if(traverse1)
        {
            return tableIterator2d(data,0,idx2,traverse1);
        }
        else{
            
            return tableIterator2d(data,idx1,0,traverse1);
        }
    }
    tableIterator2d end()
    {
        if(traverse1)
        {
            return tableIterator2d(data,data->size()-1,idx2,traverse1);
        }
        else{
            
            return tableIterator2d(data,idx1,(*data)[idx1]->size()-1,traverse1);
        }
    }
    void operator=(const tableIterator2d& other)
    {
        this->data = other.data;
        this->idx1 = other.idx1;
        this->idx2 = other.idx2;
        this->traverse1 = other.traverse1;
    }
    scalar operator*()
    {
        return (*data)[idx1][idx2];   
    }
    tableIterator2d& operator++()
    {
        idx++;
        return *this;
    }
    scalar index() {return idx;}
    friend bool operator==(tableIterator2d lhs, tableIterator2d rhs) 
    {
        return ((lhs.data) == (rhs.data)) 
                && (lhs.idx1 == rhs.idx1) 
                && (lhs.idx2 == rhs.idx2) 
                && (lhs.traverse1 == rhs.traverse1);
    }
    friend bool operator!=(tableIterator2d lhs, tableIterator2d rhs)  
    {
        return !(lhs==rhs);
    }
    
};

#endif /* PROPELLER_INTERPOLATION_REGULAR_LINEAR_ITERATORS */
