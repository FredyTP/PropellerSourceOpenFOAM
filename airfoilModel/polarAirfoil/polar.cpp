#include "polar.H"
#include "closestNeighbor.H"

foam::polar::polar(List<scalar> alpha, List<scalar> cl, List<scalar> cd, scalar Re, scalar Ma)
{
    cl_alfa = new closestNeighbor(alfa,cl);
    cd_alfa = new closestNeighbor(alfa,cd);

}