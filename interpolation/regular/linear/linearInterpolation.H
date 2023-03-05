#ifndef PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEARINTERPOLATION1
#define PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEARINTERPOLATION1

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
private:
    FixedList<List<typeIn>,dim> inputs;
    List<typeOu> outputs;
    FixedList<label,dim> size;
public:

    /**
     *  Create interpolation class from ordered inputs and outputs ->  
     *  input: ((1 2 3)  (1 2 3)) output: (1 2 3 4 5 6 7 8 9)
     * This input would generate the following table:
     *    input  |  output
     *   ---------------------
     *   (1 1)  ->   x1
     *   (2 1)  ->   x2
     *   (3 1)  ->   x3
     *   (1 2)  ->   x4
     *   (2 2)  ->   x5
     *   (3 2)  ->   x6
     *   (1 3)  ->   x7
     *   (2 3)  ->   x8
     *   (3 3)  ->   x9
     * 
     * 
    */
    linearInterpolation(FixedList<List<typeIn>,dim>> inputs_, List<typeOu> outputs_) 
        : inputs(inputs_), outputs(outputs_)
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
    linearInterpolation(List<FixedList<typeIn>,dim>> inputs_, List<typeOu> outputs_,List<dim> size) 
        : inputs(inputs_), outputs(outputs_)
    {
        if(!checkSize())
        {
            Info<<"interpolation sizes are inconsistent"<<endl;
        }
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
    label index(FixedList<label,dim> idx)
    {
        label p = 0; // i + j*N1 + k*N1*N2 + w*N1*N2*N3 ...
        for(label d = 0 ; d<dim; d++)
        {
            label i_n = idx[d]; // indice de la dimension d
            label tam = 1;
            for(label d2 = 0;d2<d;d2++)
            {
                tam*=size[d2];
            }
            p+=i_n*tam;
        }
        return p;
    }
};


}


#endif //PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEARINTERPOLATION1
