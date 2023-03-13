#ifndef PROPELLER_INTERPOLATION_REGULAR_REGULAR_INTERPOLATION
#define PROPELLER_INTERPOLATION_REGULAR_REGULAR_INTERPOLATION

#include "List.H"
#include "scalar.H"
#include "interpolationTable.H"
#include "interpolated.H"

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
class regularInterpolation : public interpolationTable<typeIn,typeOu,dim>
{

public:

    /**
     * @brief returns a list of interpolated coefficients and nodes
     * 
     * The reconstructed interpolated value is : coeff_1 * value_1 + coeff_2 * value_2 ...
    */
    virtual interpolated<typeIn,typeOu> interpolate(FixedList<typeIn,dim> input_) = 0;

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
            tam *=size[d];
        }
        return tam;
    }
    friend Ostream& operator<<(Ostream& os, const regularInterpolation<typeIn,typeOu,dim>& table)
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
protected:

    //Inputs of the N-dimensional cube iterated in specified order
    FixedList<List<typeIn>,dim> inputs;

    //Ordered outputs in iterating order
    List<typeOu> outputs;

    //Size of each dimension representing the discretizations of Ns-dimensional cube
    FixedList<label,dim> size;
    //Iteration order in which the dimensions are traveled (0 1 2 3...)
    FixedList<label,dim> order;
public:
    regularInterpolation()
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

    regularInterpolation(FixedList<List<typeIn>,dim> inputs_, List<typeOu> outputs_) 
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
    regularInterpolation(FixedList<List<typeIn>,dim> inputs_, List<typeOu> outputs_,FixedList<label,dim> iterationOrder) 
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
    /*regularInterpolation(List<FixedList<typeIn,dim>> inputs_, List<typeOu> outputs_,bool raw) 
    {
        this->setRawData(inputs_,outputs_);
    }*/
    void setRawData(List<List<typeIn>>& inputs_, List<typeOu>& outputs_)
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
                for(int c = 0; c <inputs[d].size(); ++c)
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

        Info<<"Unsorted inputs"<<endl;
        for(label d=0; d< dim; d++)
        {
            Info<<inputs[d]<<endl;
        }
        
        //Sort arrays
        Info<< "Sorted inputs"<<endl;
        for(label d=0; d< dim; d++)
        {
            std::sort(inputs[d].begin(),inputs[d].end());
            size[d]=inputs[d].size();
            Info<<inputs[d]<<endl;
        }

        outputs.resize(this->getOutputSize());
        
        //Fill outputs in place
        for(label i = 0 ; i< outputs_.size();i++)
        {
            FixedList<label,dim> idx;
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
        
        Info<<"Outputs: "<<endl;
        Info<<outputs<<endl;
        //Fill gaps with linear interpolation
    }
    bool checkSize()
    {
        for(label i = 0; i< dim; i++)
        {
            size[i]=inputs[i].size();
        }
        label len = std::accumulate(size.begin(),size.end(),1,[](label a, label b){return a*b;});
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
                tam*=size[order[d2]];
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
                tam*=size[order[d2]];
            }
            idxs[order[d]]=p/tam;

            p-= idxs[order[d]]*tam;
        }
        return idxs;
    }

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



}

#endif /* PROPELLER_INTERPOLATION_REGULAR_REGULAR_INTERPOLATION */
