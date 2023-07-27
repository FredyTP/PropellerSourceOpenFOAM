#ifndef PROPELLER_INTERPOLATION_REGULAR_REGULAR_INTERPOLATION
#define PROPELLER_INTERPOLATION_REGULAR_REGULAR_INTERPOLATION

#include "List.H"
#include "scalar.H"
#include "InterpolationTable.H"
#include "Interpolated.H"

namespace Foam
{

/**
 * @brief base class for regular interpolation
 * Implements base index 1D -> n-D mapping and n-D -> 1D mapping
 * 
 * It is possible to make a N-dimesional regular interpolation
 * using a 1D array with a N-dimesional index mapping
 * index = N*i + j;
 * index = Nx*Ny*i + Nx*j + k; and so on...
 * 
 * the number of coefficients can be "easily" obtained
*/
template<class typeIn, class typeOu, label dim>
class RegularInterpolation : public InterpolationTable<typeIn,typeOu,dim>
{

public:

    virtual ~RegularInterpolation() = default;
    
    /**
     * @brief returns a list of Interpolated coefficients and nodes
     * 
     * The reconstructed Interpolated value is : coeff_1 * value_1 + coeff_2 * value_2 ...
    */
    virtual Interpolated<typeIn,typeOu> interpolate(const FixedList<typeIn, dim>& input_) const = 0;

    typeOu getOutput(label idx)
    {
        return outputs[idx];
    }
    typeOu getOutput(FixedList<label,dim> idxs)
    {
        return outputs[index(idxs)];
    }
    const typeOu getOutput(label idx) const
    {
        return outputs[idx];
    }
    const typeOu getOutput(FixedList<label,dim> idxs) const
    {
        return outputs[index(idxs)];
    }
    scalar getOutputSize() const
    {
        scalar tam = 1;
        for(label d =0; d<dim;d++)
        {
            tam *=size_[d];
        }
        return tam;
    }
    friend Ostream& operator<<(Ostream& os, const RegularInterpolation<typeIn,typeOu,dim>& table) 
    {
        for(label i=0;i<table.outputs.size();i++)
        {
            FixedList<label,dim> idx = table.invIndex(i);
            os<<"(";
            for(label j=0;j<dim;j++)
            {
                if(j!=0)
                {
                    os<<", ";
                }
                os<<table.inputs[j][idx[j]];
            }
            os<<") -> " <<table.outputs[i]<<endl;
        }
        return os;
    }

    const FixedList<List<typeIn>,dim>& getInputs() const { return inputs; }
    const List<typeOu>& getOutputs() const { return outputs;}
protected:

    //Inputs of the N-dimensional cube iterated in specified order
    FixedList<List<typeIn>,dim> inputs;

    //Ordered outputs in iterating order
    List<typeOu> outputs;

    //Size of each dimension representing the discretizations of Ns-dimensional cube
    FixedList<label,dim> size_;
    //Iteration order in which the dimensions are traveled (0 1 2 3...)
    FixedList<label,dim> order;
public:
    RegularInterpolation()
    {
        
    }

    void setData(FixedList<List<typeIn>,dim> inputs_, List<typeOu> outputs_)
    {
        inputs = inputs_;
        outputs = outputs_;
        for(label i = 0; i < dim; i++)
        {
            order[i]=i;
        }
        if(!checkSize())
        {
            Info<<"interpolation sizes are inconsistent"<<endl;
        }
    }
    /**
     *  Create interpolation class from ordered inputs and outputs ->  
     *  input: ((1 2 3)  (1 2 3)) output: (1 2 3 4 5 6 7 8 9)
     * This input would generate the following table: 
     * First table with iterating ording by default and second with order (1 0)
     * 
     * 
     *     default (0 1)          order (1 0)
     *    input  |  output      input  |  output
     *   -----------------    --------------------
     *   (1 1)  ->   x1         (1 1)  ->   x1
     *   (2 1)  ->   x2         (1 2)  ->   x2
     *   (3 1)  ->   x3         (1 3)  ->   x3
     *   (1 2)  ->   x4         (2 1)  ->   x4
     *   (2 2)  ->   x5         (2 2)  ->   x5
     *   (3 2)  ->   x6         (2 3)  ->   x6
     *   (1 3)  ->   x7         (3 1)  ->   x7
     *   (2 3)  ->   x8         (3 2)  ->   x8
     *   (3 3)  ->   x9         (3 3)  ->   x9
     *  
    */

    RegularInterpolation(FixedList<List<typeIn>,dim> inputs_, List<typeOu> outputs_) 
        : inputs(inputs_), outputs(outputs_)
    {
        for(label i = 0; i < dim; i++)
        {
            order[i]=i;
        }
        if(!checkSize())
        {
            Info<<"interpolation sizes are inconsistent"<<endl;
        }
    }
    RegularInterpolation(FixedList<List<typeIn>,dim> inputs_, List<typeOu> outputs_,FixedList<label,dim> iterationOrder) 
        : inputs(inputs_), outputs(outputs_), order(iterationOrder)
    {
        
        if(!checkSize())
        {
            Info<<"interpolation sizes are inconsistent"<<endl;
        }
    }
    /**
     *                                                                   input  |  output
     *                                                                ---------------------
     * Create interpolation class from raw inputs points and outputs   (1 1 1)  ->   x1
     *                                                                 (2 1 1)  ->   x2
     *                      size: (2 2 2)                              (1 2 1)  ->   x3
     *                                                                 (2 2 1)  ->   x4
     *                                                                 (1 1 2)  ->   x5
     *                                                                 (2 1 2)  ->   x6
     *                                                                 (1 2 2)  ->   x7
     *                                                                 (2 2 2)  ->   x8
    */
    virtual bool setRawData(const List<List<typeIn>>& inputs_, const List<typeOu>& outputs_)
    {
        for(label i = 0; i < dim; i++)
        {
            order[i]=i;
        }

        //Find all diferent variables
        for(label i = 0; i < inputs_.size(); i++)
        {
            for(label d=0; d< dim; d++)
            {
                bool contains = false;
                //Check if that input already been processed
                for(int c = 0; c <inputs[d].size(); ++c) //Refactor with list.find()
                {
                    if(inputs[d][c]==inputs_[i][d])
                    {
                        contains=true;
                        break;
                    }
                }
                if(!contains)
                {
                    inputs[d].push_back(inputs_[i][d]);
                }
            }
        }

        //Sort arrays

        for(label d=0; d< dim; d++)
        {
            std::sort(inputs[d].begin(),inputs[d].end());
            size_[d]=inputs[d].size();
        }

        //Output size != N1 * N2 * ...Ni; not regular gridded data
        if(outputs_.size() != this->getOutputSize())
        {
            //This criteria is not 100% robust but easy check
            return false;
        }
        
        outputs.resize(this->getOutputSize());
        
        //Fill outputs in place
        for(label i = 0 ; i< outputs_.size();i++)
        {
            FixedList<label,dim> idx(0);
            for(label d = 0; d<dim;d++)
            {              
                for(label j = 0; j < inputs[d].size(); j++)
                {
                    if(inputs[d][j]==inputs_[i][d])
                    {
                        idx[d]=j;
                    }
                }
            }
            outputs[this->index(idx)]=outputs_[i];
        }
        
        return true;
        //Fill gaps with linear interpolation
    }
    bool checkSize()
    {
        for(label i = 0; i< dim; i++)
        {
            size_[i]=inputs[i].size();
        }
        label len = std::accumulate(size_.begin(),size_.end(),1,[](label a, label b){return a*b;});
        label oLen = outputs.size();
        return len == oLen;
    }
    label index(FixedList<label,dim> idx) const
    {
        label p = 0; // i + j*N1 + k*N1*N2 + w*N1*N2*N3 ...
        for(label d = 0 ; d<dim; d++)
        {
            label i_n = idx[order[d]]; // indice de la dimension d
            label tam = 1;
            for(label d2 = 0;d2<d;d2++)
            {
                tam*=size_[order[d2]];
            }
            p+=i_n*tam;
        }
        return p;
    }
    FixedList<label,dim> invIndex(label index) const 
    {
        label p = index; // i + j*N1 + k*N1*N2 + w*N1*N2*N3 ...
        FixedList<label,dim> idxs;
        
        for(label d = dim-1 ; d>=0; --d)
        {
            label tam = 1;
            for(label d2 = 0;d2<d;d2++)
            {
                tam*=size_[order[d2]];
            }
            idxs[order[d]]=p/tam;

            p-= idxs[order[d]]*tam;
        }
        return idxs;
    }

    /**
     * @brief Find linear interpolation indexes
     * 
     * @param input value that will be Interpolated
     * @param inputs gridded data points
     * @param i1 index of closest data point < input
     * @param i2 index of closest data point > input
     * @return result, if result = 0, input is out of bounds
     * if result = 1 data is inside bounds
     * 
     * if i1 = i2, input value == inputs[i1] (i1 and i2)
     */
    static label FindIndex(typeIn input, const List<typeIn>& inputs, label &i1, label &i2)
    {
        if(input<inputs[0])
        {
            //indexes are both set to closest index
            i1=0;
            i2=0;
            //out of bounds
            return 0;
        }
        for(label i = 0; i< inputs.size()-1; i++)
        {
            //If (luckily) there is coincidence
            if(input == inputs[i])
            {
                i1=i;
                i2=i;
                return 1;
            }

            //Between i and i+1
            if(inputs[i] < input && input< inputs[i+1])
            {
                i1=i;
                i2=i+1;
                return 1;
            }
        }

        if(input==inputs[inputs.size()-1])
        {
            i1=inputs.size()-1;
            i2=i1;
            return 1;
        }

        //Outside upper bounds
        i1=inputs.size()-1;
        i2=i1;
        return 0;

    }

};

template<class typeIn, class typeOu, label dim>
class ConstantInterpolation : public RegularInterpolation<typeIn,typeOu,dim>
{
    public:
    ConstantInterpolation(const typeOu& value) : value_(value) {}

    bool setRawData(List<List<typeIn>>& inputs_, List<typeOu>& outputs_) { return true;}
    Interpolated<typeIn,typeOu> interpolate(const FixedList<typeIn, dim>& input) const
    {
        return Interpolated<typeIn,typeOu>({1},{value_});
    }
    label size() { return 1; }
    private:
    typeOu value_;
};

}

#endif /* PROPELLER_INTERPOLATION_REGULAR_REGULAR_INTERPOLATION */
