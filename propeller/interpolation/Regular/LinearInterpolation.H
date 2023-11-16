#ifndef PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEARINTERPOLATION
#define PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEARINTERPOLATION

#include "List.H"
#include "scalar.H"
#include "RegularInterpolation.H"


namespace Foam
{


/**
 * @brief Interpolation class for N-dimensional linear interpolation
 * 
 * @tparam typeIn 
 * @tparam typeOu 
 * @tparam dim 
 */
template<class typeIn, class typeOu, label dim>
class LinearInterpolation : public RegularInterpolation<typeIn,typeOu,dim>
{
private:
    typedef RegularInterpolation<typeIn,typeOu,dim> parent;
public:

    LinearInterpolation()
    {
        
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

    LinearInterpolation(FixedList<List<typeIn>,dim> inputs_, List<typeOu> outputs_) 
        : RegularInterpolation<typeIn,typeOu,dim>(inputs_,outputs_)
    {

    }
    LinearInterpolation(FixedList<List<typeIn>,dim> inputs_, List<typeOu> outputs_,FixedList<label,dim> iterationOrder) 
        : RegularInterpolation<typeIn,typeOu,dim>(inputs_,outputs_,iterationOrder)
    {

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
    LinearInterpolation(List<FixedList<typeIn,dim>> inputs_, List<typeOu> outputs_,bool raw) 
        : RegularInterpolation<typeIn,typeOu,dim>(inputs_,outputs_,raw)
    {

    }
    
    Interpolated<typeIn,typeOu> interpolate(const FixedList<typeIn, dim>& input) const
    {
        Interpolated<typeIn,typeOu> iReturn;
        auto& coeff = iReturn.coefficients();
        auto& values = iReturn.points();


        FixedList<FixedList<label,2>,dim> indexIt;
        FixedList<label,dim> result;
        FixedList<typeIn,dim> delta;

        //Get index and results
        for(label i = 0; i< dim; i++)
        {
            result[i] = RegularInterpolation<typeIn,typeOu,dim>::FindIndex(input[i],parent::inputs[i],indexIt[i][0],indexIt[i][1]);
            //<out of bounds>
            if(result[i]==0)
            {                
                switch (this->exMode)
                {
                case extrapolationMode::emConstant :
                    delta[i] = 0.5;
                    break;
                case extrapolationMode::emZero :
                    delta[i] = 0;
                    break;
                default:
                    delta[i] = 0;
                    break;
                }
                
            }
            //<coincident>
            else if(indexIt[i][0]==indexIt[i][1])
            {
                delta[i] = 0.5;
            }
            else{
                delta[i] = (input[i] - parent::inputs[i][indexIt[i][0]])/(parent::inputs[i][indexIt[i][1]]-parent::inputs[i][indexIt[i][0]]);
            }
        }

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
                if(result[j]==0)
                {
                    coeff[i] *= delta[j];
                }
                else
                {
                    coeff[i]*=(1-z[j]+sign*delta[j]);
                }
            };
            values[i]= parent::outputs[parent::index(idx)];

            //Iterate all permutations in binary secuence
            for (int k = 0; k < dim; ++k) 
            {
                if (z[k] != 1) {
                    ++z[k];
                    break;
                } else {
                    z[k] = 0;
                }
            }
        }

        return iReturn;
    }    

    virtual label size() {
        label s = 1;
        forAll(this->size_,i)
        {
            s*=this->size_[i];
        }
        return s;
    }
};


}

#endif //PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEARINTERPOLATION

