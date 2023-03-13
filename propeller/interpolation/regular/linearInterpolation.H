#ifndef PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEARINTERPOLATION
#define PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEARINTERPOLATION

#include "List.H"
#include "scalar.H"
#include "regularInterpolation.H"


namespace Foam
{


template<class typeIn, class typeOu, label dim>
class linearInterpolation : public regularInterpolation<typeIn,typeOu,dim>
{
private:
    typedef regularInterpolation<typeIn,typeOu,dim> parent;
public:

    linearInterpolation()
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

    linearInterpolation(FixedList<List<typeIn>,dim> inputs_, List<typeOu> outputs_) 
        : regularInterpolation<typeIn,typeOu,dim>(inputs_,outputs_)
    {

    }
    linearInterpolation(FixedList<List<typeIn>,dim> inputs_, List<typeOu> outputs_,FixedList<label,dim> iterationOrder) 
        : regularInterpolation<typeIn,typeOu,dim>(inputs_,outputs_,iterationOrder)
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
    linearInterpolation(List<FixedList<typeIn,dim>> inputs_, List<typeOu> outputs_,bool raw) 
        : regularInterpolation<typeIn,typeOu,dim>(inputs_,outputs_,raw)
    {

    }
    
    interpolated<typeIn,typeOu> interpolate(FixedList<typeIn,dim> input) const
    {
        interpolated<typeIn,typeOu> iReturn;
        auto& coeff = iReturn.coefficients();
        auto& values = iReturn.points();


        FixedList<FixedList<label,2>,dim> indexIt;
        FixedList<label,dim> result;
        FixedList<typeIn,dim> delta;

        //Get index and results
        for(label i = 0; i< dim; i++)
        {
            result[i] = regularInterpolation<typeIn,typeOu,dim>::FindIndex(input[i],parent::inputs[i],indexIt[i][0],indexIt[i][1]);
            if(result[i]==0 || indexIt[i][0]==indexIt[i][1])
            {
                //<out of bounds> || <coincident>
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
                coeff[i]*=(1-z[j]+sign*delta[j]);
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
};


}

#endif //PROPELLER_INTERPOLATION_REGULAR_LINEAR_LINEARINTERPOLATION

