#ifndef PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEARINTERPOLATION
#define PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEARINTERPOLATION

#include "List.H"
#include "scalar.H"
#include "regularInterpolation.H"
namespace Foam
{

class linearInterpolation1D 
: public regularInterpolation1D
{
protected:
    const List<scalar> inputs_;
    const List<scalar> outputs_;
public:

    //Constructors

        //- Construct from input and output list
        linearInterpolation1D
        (
            const List<scalar> inputs,
            const List<scalar> outputs
        );
        
    virtual scalar interpolate(FixedList<scalar,1> input);

    
protected:
    /**
     * @return 0 if out of bounds 1 if succesful
    */
    label findIndex(scalar input,label& i1, label& i2) const;

    void interpolationCoefficients(scalar input, label i1, label i2, scalar& a1, scalar &a2) const;
    
public:
   
    static void InterpolationCoefficients(scalar input,const List<scalar>& inputs, label i1, label i2, scalar& a1, scalar &a2);
    static label FindIndex(scalar input, const List<scalar>& inputs, label &i1, label &i2);
    
};

template<class typeIn, class typeOu, label dim>
class linearInterpolation
{
public:
    typeOu getOutput(label idx)
    {
        return outputs[idx];
    }
    typeOu getOutput(FixedList<label,dim> idxs)
    {
        return outputs[index(idxs)];
    }
private:

    //Inputs of the N-dimensional cube iterated in specified order
    FixedList<List<typeIn>,dim> inputs;

    //Ordered outputs in iterating order
    List<typeOu> outputs;

    //Size of each dimension representing the discretizations of N-dimensional cube
    FixedList<label,dim> size;
    //Iteration order in which the dimensions are traveled (0 1 2 3...)
    FixedList<label,dim> order;
public:

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
    linearInterpolation(FixedList<List<typeIn>,dim> inputs_, List<typeOu> outputs_) 
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
    linearInterpolation(FixedList<List<typeIn>,dim> inputs_, List<typeOu> outputs_,FixedList<label,dim> iterationOrder) 
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
    /*linearInterpolation(List<FixedList<typeIn,dim>> inputs_, List<typeOu> outputs_,FixedList<label,dim> size) 
        : inputs(inputs_), outputs(outputs_), order({0,1,2})
    {
        if(!checkSize())
        {
            Info<<"interpolation sizes are inconsistent"<<endl;
        }
    }*/
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
    label index(FixedList<label,dim> idx)
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

    void interpolationList(FixedList<typeIn,dim> input,List<typeIn>& coeff,List<typeOu>& values)
    {
        FixedList<FixedList<label,2>,dim> indexIt;
        FixedList<label,2> result;
        FixedList<typeIn,2> delta;

        //Get index and results
        for(label i = 0; i< dim; i++)
        {
            result[i] = linearInterpolation1D::FindIndex(input[i],inputs[i],indexIt[i][0],indexIt[i][1]);
            if(result[i]==0 || indexIt[i][0]==indexIt[i][1])
            {
                //<out of bounds> || <coincident>
                delta[i] = 0.5;
            }
            else{
                delta[i] = (input[i] -inputs[i][indexIt[i][0]])/(inputs[i][indexIt[i][1]]-inputs[i][indexIt[i][0]]);

            }
        }
        Info<<"deltas" <<delta<<endl;

        //Resize to fit all coefficients
        label numrow = pow(2,dim);
        coeff.resize(numrow);
        values.resize(numrow);

        std::vector<int> z(dim, 0);
        //Init coeffs
        for(label i = 0; i< numrow ;i++)
        {
            coeff[i]=1;
        }

        //Init indexes
        for(label i=0;i<numrow;i++)
        {

            FixedList<label,dim> idx;
            for(label j=0;j<dim;j++)
            {
                idx[j]=indexIt[j][z[j]];               
                scalar sign = z[j]==0? -1: 1;
                coeff[i]*=(1-z[j]+sign*delta[j]);
            }
            Info<<idx<<endl;
            values[i]=outputs[index(idx)];

            //Iterate all permutations in binary secuence
            for (int k = 0; k < 2; ++k) 
            {
                if (z[k] != 1) {
                    ++z[k];
                    break;
                } else {
                    z[k] = 0;
                }
            }
        }

    }
};


}


#endif //PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEARINTERPOLATION
