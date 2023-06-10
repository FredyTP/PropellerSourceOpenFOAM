#ifndef PROPELLER_UTIL_FUNCTION_SOLVER
#define PROPELLER_UTIL_FUNCTION_SOLVER

#include <functional>
#include "scalar.H"
#include "FixedList.H"
#include "simpleMatrix.H"

namespace Foam 
{
namespace util
{


/**
 * @brief Wrapper template class for function 0 finder.
 * The available type of solver is for N-dim input N-dim output functions
 *  
 */
class functionSolver
{

public:
    typedef scalarField varType;
    typedef std::function<varType(varType)> funcType;


    static varType NewtonRapson
    (
        label N,
        const funcType& f,
        const varType& x0,
        const varType& dx,
        scalar relax,
        label maxIter,
        scalar tol,
        bool verbose
    );


};

inline typename functionSolver::varType functionSolver::NewtonRapson
(
    label N,
    const funcType& f,
    const varType& x0,
    const varType& dx,
    scalar relax,
    label maxIter,
    scalar tol,
    bool verbose
)
{
    scalar err = GREAT;
    label iter = 0;

    SquareMatrix<scalar> J(N,Zero);

    varType old(N,Zero);

    varType x(x0);

    while((err>tol) && (iter < maxIter))
    {
        varType xlast(x);
        old = f(x);

        for(label i = 0;i<N;i++)
        {
            x[i] -= dx[i]/2.0;
            varType y0 = f(x);

            x[i] += dx[i];
            varType y1 = f(x);

            varType dy = (y1-y0)/dx[i];

            //Set J cmp
            for(label j =0; j <N;j++)
            {
                J(j,i) = dy[j];
            }
            x=xlast;
        }
        solve(J,old);
        
        varType xNew = x - relax * old;

        varType ytest = f(xNew);
        err = sqrt(sumSqr(ytest));

        x = xNew;
        ++iter;
    }

    if(verbose)
    {
        if (iter == maxIter)
        {
            Info<< "    solution not converged in " << iter
                << " iterations, final residual = " << err
                << "(" << tol << ")" << endl;
        }
        else
        {
            Info<< "    final residual = " << err << "(" << tol
                << "), iterations = " << iter << endl;
        }
    }

    return x;

}


}
}
#endif /* PROPELLER_UTIL_FUNCTION_SOLVER */
