/*==========================================================
 *
 * align_fast.cpp
 *
 * Aligns two sequences of data and returns the resulting optimal path
 *
 *========================================================*/

#include "matrix.h"
#include <stdlib.h>
#include "mex.h"
#include <math.h>

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mxArray* res = 0;
    double d = 0;
    int i;
    
    mexPrintf("Is struct? %d\n",mxIsStruct(prhs[0]));
    
    res = mxGetField(prhs[0],0,"level_mean");
    if (res == NULL)
    {
        mexPrintf("No field level_mean found!\n");
    }
    else
    {
        mexPrintf("Size of level_mean: %d\n",mxGetNumberOfElements(res));
    }
    
    for (i=0; i<1e6; i++)
    {
        d += log(-5 + (2&3));
    }
    mexPrintf("Foo %f\n",d);
}
