/*==========================================================
 *
 * viterbi1d.cpp
 *
 * Runs Viterbi algorithm on a single level sequence
 *
 * [dpath] = viterbi1d(event);
 *
 *========================================================*/

#include "matrix.h"
#include "mex.h"
#include <stdint.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "viterbi.h"

using namespace std;

// skip and stay probabilities
const double skip_prob = 0.10;
const double stay_prob = 0.01;

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

    EventData event = EventData(prhs[0]);
    
    int n = event.length;
    
    // create a transition matrix
    double* T = new double[N_STATES*N_STATES];
    buildT(T, skip_prob, stay_prob);
    
    // and logify it
    for (int i=0; i<N_STATES*N_STATES; i++)
        T[i] = log(1/inf + T[i]);
    
    // allocate n x N_STATES array to hold likelihoods
    // and another to hold back-pointers
    double *liks = new double[n*N_STATES];
    int *backptrs = new int[n*N_STATES];

    mexPrintf("Initialized, starting iterations.\n");
    
    // seed initial level with observation likelihoods, no transitions
    // could normalize, but who cares?
    for (int i=0; i<N_STATES; i++)
    {
        // use per-state deviation
        liks[i] = log(1/inf + event.model.getProb(event.mean[0],event.stdv[0],i));
        backptrs[i] = -1;
    }
    
    // now loop through each level
    for (int i=1; i<n; i++)
    {
        // this observation's likelihood bit
        double* curlik = liks + i*N_STATES;
        
        // calculate log-normalized observations
        double obs[N_STATES];
        
        for (int j=0; j<N_STATES; j++)
            obs[j] = log(1/inf + event.model.getProb(event.mean[i],event.stdv[i],j));
        
        for (int curst=0; curst<N_STATES; curst++)
        {
            // this column of transition matrix
            double* Tcol = T+curst*N_STATES;

            double lmax = -inf;
            int lstate = -1;
            
            for (int prevst=0; prevst<N_STATES; prevst++)
            {
                double l = obs[curst] + Tcol[prevst];
                l += curlik[prevst-N_STATES];
                if (l>lmax)
                {
                    lmax = l;
                    lstate = prevst;
                }
            }
            // now finally check stays
            double l = obs[curst]+log(stay_prob)+curlik[curst-N_STATES];
            if (l>lmax)
            {
                lmax = l;
                lstate = curst;
            }
            
            // and save likelihood
            curlik[curst] = lmax;
            backptrs[curst + i*N_STATES] = lstate;            
        }
        if (i%100 == 0)
        {
            mexPrintf("Iteration %d\n",i);
            cout.flush();
        }
    }
    
    // now output the resulty thingy
    plhs[0] = mxCreateDoubleMatrix(n,2,mxREAL);
    double *pr = mxGetPr(plhs[0]);
    
    // find largest final likelihood
    double curlik = -inf;
    int curst = -1;
    
    for (int i=0; i<N_STATES; i++)
    {
        if (liks[i+(n-1)*N_STATES] > curlik)
        {
            curlik = liks[i+(n-1)*N_STATES];
            curst = i;
        }
    }
    
    for (int i=n-1; i>=0; i--)
    {
        curlik = liks[i*N_STATES+curst];
        pr[i] = curlik;
        // 1-based states for matlab
        pr[n+i] = curst + 1;
        curst = backptrs[i*N_STATES+curst];
    }
    
    delete[] T;
    delete[] liks;
    delete[] backptrs;
}
