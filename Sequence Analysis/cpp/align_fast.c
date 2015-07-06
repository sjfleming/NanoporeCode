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
    double *data0;                  /* 1xM input matrix */
    double *data1;                  /* 1xN input matrix */
    double *dpath;
    double sd=2.0;                  /* standard dev */
    double si=0.4;                  /* insertion penalty */
    
    int n0 = mxGetNumberOfElements(prhs[0]);
    int n1 = mxGetNumberOfElements(prhs[1]);
    
    // now allocate score matrix and steps matrix
    double* scores = (double*)mxCalloc((n0+1)*(n1+1),sizeof(double));
    uint8_T* steps = (uint8_T*)mxCalloc((n0+1)*(n1+1),sizeof(uint8_T));
    double* dptemp = (double*)mxCalloc((n0+n1)*2,sizeof(double));
    
    int w = n0+1;
    
    double s[3];
    
    int i,j,k;
    
    mexPrintf("Got arrays of size %d and %d\n",n0,n1);
    
    // check for proper number of arguments
    if(nrhs<2) {
        mexErrMsgIdAndTxt("CrampBio:align_fast","At least two inputs required.");
    }
    if(nlhs>2) {
        mexErrMsgIdAndTxt("CrampBio:align_fast","At most two outputs.");
    }
    
    // make sure the second input argument is type double
    for (i=0; i<2; i++)
    {
        if( !mxIsDouble(prhs[i]) || 
             mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt("CrampBio:align_fast","Input matrices must be type double.");
        }
    }
    
    if (nrhs >= 3)
        sd = mxGetScalar(prhs[2]);
    if (nrhs >= 4)
        si = mxGetScalar(prhs[3]);
    
    
    data0 = mxGetPr(prhs[0]);
    data1 = mxGetPr(prhs[1]);

    scores[0] = 0.0;
    steps[0] = 3;
    
    mexPrintf("Initialized, starting iterations.\n");
    
    
    // 1 is from above, 2 is from the left, 3 is diag
    for (i=1; i<=n0; i++)
    {
        scores[i] = i*si;
        steps[i] = 1;
    }
    for (j=1; j<=n1; j++)
    {
        scores[j*w] = j*si;
        steps[j*w] = 2;
    }
    
    mexPrintf("Starting scoring.\n");
    
    // fill up score matrix
    for (j=0; j<n1; j++)
    {
        for (i=0; i<n0; i++)
        {
            int mi = 1;
            // find minimal score out of the three up and to the left
            double d = (data0[i]-data1[j])*1e3/sd;
            // up and to the left, match score
            s[2] = (1-exp(-0.5*d*d)) + scores[i + w*j];
            // to the left, insertion score
            s[0] = si + scores[i + w*(j+1)];
            // up, insertion score
            s[1] = si + scores[(i+1) + w*j];
            
            if (s[1] < s[0])
            {
                s[0] = s[1];
                mi = 2;
            }
            if (s[2] < s[0])
            {
                s[0] = s[2];
                mi = 3;
            }
            
            scores[(i+1) + w*(j+1)] = s[0];
            steps[(i+1) + w*(j+1)] = mi;
        }
    }
    
    mexPrintf("Scoring complete.\n");

    // Now reconstruct the sequence that led here, as indices

    i = n0+1;
    j = n1+1;
    k = 0;
    
    while (i>0 && j>0)
    {
        int st = steps[(i-1) + w*(j-1)];
        
        dptemp[2*k] = i;
        dptemp[2*k+1] = j;
        
        switch (st)
        {
            case 1:
                i = i - 1;
                break;
            case 2:
                j = j - 1;
                break;
            case 3:
                i = i - 1;
                j = j - 1;
                break;
        }
        k = k + 1;
    }

    mexPrintf("Aligned with %d total steps\n",k);
    
    // create output dpath
    plhs[0] = mxCreateDoubleMatrix(k,2,mxREAL);
    dpath = mxGetPr(plhs[0]);
    
    for (i=0; i<k; i++)
    {
        dpath[i] = dptemp[2*(k-i-1)];
        dpath[i+k] = dptemp[2*(k-i-1)+1];
    }
    
    // create the output for scores array too
    if (nlhs > 1)
    {
        // create size-zero matrix to avoid allocating memory
        plhs[1] = mxCreateDoubleMatrix(0,0,mxREAL);
        // set pointer and size
        mxSetPr(plhs[1], scores);
        mxSetM(plhs[1],n0+1);
        mxSetN(plhs[1],n1+1);
        // it will now get returned
    }
    else
    {
        // otherwise, delete it
        mxFree(scores);
    }
    
    mxFree(steps);
    mxFree(dptemp);
}
