#ifndef PROPELLER_INTERPOLATION_REGULAR_QUADRATICINTERPOLATION
#define PROPELLER_INTERPOLATION_REGULAR_QUADRATICINTERPOLATION

#include "List.H"
#include "scalar.H"
#include "regularInterpolation.H"


namespace Foam
{

/**
 * @brief Regular quadratic order interpolation class
 * @remark Interpolation data should be of size >= 3 in each dimension,
 * if the data set is smaller, consider linear interpolation
 * 
 * @tparam typeIn 
 * @tparam typeOu 
 * @tparam dim 
 */
template<class typeIn, class typeOu, label dim>
class quadraticInterpolation : public regularInterpolation<typeIn,typeOu,dim>
{
private:
    typedef regularInterpolation<typeIn,typeOu,dim> parent;
public:

    quadraticInterpolation()
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

    quadraticInterpolation(FixedList<List<typeIn>,dim> inputs_, List<typeOu> outputs_) 
        : regularInterpolation<typeIn,typeOu,dim>(inputs_,outputs_)
    {

    }
    quadraticInterpolation(FixedList<List<typeIn>,dim> inputs_, List<typeOu> outputs_,FixedList<label,dim> iterationOrder) 
        : regularInterpolation<typeIn,typeOu,dim>(inputs_,outputs_,iterationOrder)
    {

    }

    /** TODO: CHANGE THIS COMMENTS
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
    quadraticInterpolation(List<FixedList<typeIn,dim>> inputs_, List<typeOu> outputs_,bool raw) 
        : regularInterpolation<typeIn,typeOu,dim>(inputs_,outputs_,raw)
    {

    }
    

    static FixedList<typeIn,3> LagrangePolinomial(const FixedList<typeIn,3>& xi,typeIn x)
    {
        FixedList<typeIn,3> f;

        f[0] = (x-xi[1])*(x-xi[2])/(xi[0]-xi[1])/(xi[0]-xi[2]);
        f[1] = (x-xi[0])*(x-xi[2])/(xi[1]-xi[0])/(xi[1]-xi[2]);
        f[2] = (x-xi[0])*(x-xi[1])/(xi[2]-xi[0])/(xi[2]-xi[1]);

        return f;
    }
    /**
     * @brief Find quadratic interpolation indexes
     * 
     * @param input value that will be interpolated
     * @param inputs gridded data points
     * @param index sorted index in ascending input value order
     * @return result, if result = 0, input is out of bounds
     * if result = 1 data is inside bounds
     * 
     * if inside bounds
     * if index[0] = index[1], input value == inputs[i1] (i1 and i2)
     * 
     * else colapsed to closest point
     */
    static label FindIndexQuadratic(const typeIn& input, const List<typeIn>& inputs, FixedList<label,3>& index)
    {
        label i1,i2;
        label result = parent::FindIndex(input,inputs,i1,i2);
        if(i1 == i2) //either coincident or outside bounds
        {
            index[0]=i1;
            index[1]=i1;
            index[2]=i1;

            return result; //either 1 or 0
        }
        else
        {
            typeIn dist0,dist3;
            label i0,i3;
            i0 = i1-1;
            i3 = i2+1;
            if(i1 == 0) //i1 is the first point in set
            {
                //This can be done because inputs size >=3
                index[0]=i1;
                index[1]=i2;
                index[2]=i3;
                return 1; //OK
            }
            if(i2==inputs.size()-1) //i2 is the last point in set
            {
                //This can be done because inputs size >=3
                index[0]=i0;
                index[1]=i1;
                index[2]=i2;
                return 1; //OK
            }
            //Since i1 and i2 arent extreme points, data size >=4 so this is valid
            dist0 = input-inputs[i0];
            dist3 = inputs[i3]-input;
            if(dist0 < dist3) //Input is closer to i0
            {
                index[0]=i0;
                index[1]=i1;
                index[2]=i2;
                return 1; //OK
            }
            else //Input is closer to i3
            {
                index[0]=i1;
                index[1]=i2;
                index[2]=i3;
                return 1; //OK
            }
        }
    }
    interpolated<typeIn,typeOu> interpolate(FixedList<typeIn,dim> input) const
    {
        interpolated<typeIn,typeOu> iReturn;
        auto& coeff = iReturn.coefficients();
        auto& values = iReturn.points();


        FixedList<FixedList<label,3>,dim> indexIt;
        FixedList<label,dim> result;
        FixedList<FixedList<typeIn,3>,dim> f;

        //Get index and results
        for(label i = 0; i< dim; i++)
        {
            result[i] = FindIndexQuadratic(input[i],parent::inputs[i],indexIt[i]);
            //<out of bounds>
            if(result[i]==0)
            {                
                switch (this->exMode)
                {
                case extrapolationMode::emConstant :
                    f[i][0] = 1.0/3.0;
                    f[i][1] = 1.0/3.0;
                    f[i][2] = 1.0/3.0;
                    break;
                case extrapolationMode::emZero :
                    f[i][0] = 0.0;
                    f[i][1] = 0.0;                   
                    f[i][2] = 0.0;
                    break;
                default:
                    break;
                }
                
            }
            //<coincident>
            else if(indexIt[i][0]==indexIt[i][1])
            {
                f[i][0] = 1.0/3.0;
                f[i][1] = 1.0/3.0;
                f[i][2] = 1.0/3.0;
            }
            else{
                label i0,i1,i2;
                i0 = indexIt[i][0];
                i1 = indexIt[i][1];
                i2 = indexIt[i][2];

                FixedList<typeIn,3> inputValInterp={parent::inputs[i][i0],parent::inputs[i][i1],parent::inputs[i][i2]};
                f[i] = LagrangePolinomial(inputValInterp,input[i]);
            }
        }

        //Resize to fit all coefficients
        label numrow = pow(3,dim);
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
                coeff[i] *= f[j][z[j]];      
            };
            values[i]= parent::outputs[parent::index(idx)];

            //Iterate all permutations in binary secuence
            for (int k = 0; k < dim; ++k) 
            {
                if (z[k] != 2) {
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
        label s = 0;
        forAll(this->size_,i)
        {
            s+=this->size_[i];
        }
        return s;
    }
};


}

#endif //PROPELLER_INTERPOLATION_REGULAR_QUADRATICINTERPOLATION
