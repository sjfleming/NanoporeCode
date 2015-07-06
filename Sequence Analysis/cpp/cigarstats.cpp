/*==========================================================
 *
 * cigarstats.cpp
 *
 * Get statistics from a CIGAR string, reference, sequence, and start index
 * Returns total, match, mismatch, insert, deletion, clip left, clip right
 *
 *========================================================*/

#include "matrix.h"
#include "mex.h"
#include <stdint.h>
#include <cmath>
#include <iostream>

using namespace std;

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    if (nrhs != 4)
        mexErrMsgIdAndTxt("CrampBio:cigarstats","Four inputs required.");
    for (int i=0; i<3; i++)
        if (!mxIsChar(prhs[i]))
            mexErrMsgIdAndTxt("CrampBio:cigarstats","3 inputs must be type char!");
    
    uint16_t* cigar = (uint16_t*)mxGetPr(prhs[0]);
    uint16_t* ref = (uint16_t*)mxGetPr(prhs[1]);
    uint16_t* seq = (uint16_t*)mxGetPr(prhs[2]);
    
    int ncigar = mxGetNumberOfElements(prhs[0]);
    int nr = mxGetNumberOfElements(prhs[1]);
    int ns = mxGetNumberOfElements(prhs[2]);
    
    int indc = 0; // index in cigar string
    // index in reference, start at given value
    // (less one for indexing reasons)
    int indr = (int)mxGetScalar(prhs[3]) - 1;
    int inds = 0; // index in sequence
    
    int cigtot = 0; // how many total
    
    int nmatch = 0;
    int nmismatch = 0;
    int ndel = 0;
    int nins = 0;
    int nhard[2] = {0,0};
    
    // this is the current integer value of our cigar string
    int cigval = 0;
    
    // loop through cigar string
    while (indc < ncigar)
    {
        uint16_t c = cigar[indc];
        // is it a digit? if so, accumulate
        if (c >= '0' && c <= '9')
        {
            cigval = 10*cigval + (c-'0');
        }
        else
        {
            // figure out what we have
            switch (c)
            {
                case 'M':
                    // check how many matches and mismatches we have
                    for (int i=0; i<cigval; i++)
                    {
                        if (seq[inds] == ref[indr])
                            nmatch++;
                        else
                            nmismatch++;
                        
                        indr++;
                        inds++;
                    }
                    cigtot += cigval;
                    break;
                case 'I':
                    // insertion, relative to ref, so advance sequence
                    inds += cigval;
                    cigtot += cigval;
                    nins += cigval;
                    break;
                case 'D':
                    // deletion, relative to ref, so advance ref
                    indr += cigval;
                    cigtot += cigval;
                    ndel += cigval;
                    break;
                case 'S':
                    cerr << "S operation not implemented!" << endl;
                    break;
                case 'H':
                    // hard clip, keep track of this, but otherwise do nothing
                    // (just that stuff dropped from the SEQ in the bam file)
                    // is it the first hard clip?
                    if (inds==0)
                        nhard[0] = cigval;
                    else
                        nhard[1] = cigval;
                    break;
            }
            cigval = 0;
            if (inds > ns || indr > nr)
                cerr << "Index exceeded error!" << endl;
        }
        indc++;
    }
    plhs[0] = mxCreateDoubleMatrix(1,7,mxREAL);
    double *pr = mxGetPr(plhs[0]);
    pr[0] = cigtot;
    pr[1] = nmatch;
    pr[2] = nmismatch;
    pr[3] = nins;
    pr[4] = ndel;
    pr[5] = nhard[0];
    pr[6] = nhard[1];
    //cout << "Total: " << cigtot << " | Match: " << nmatch << " | Mismatch: " << nmismatch << " | Insert: " << nins << " | Del: " << ndel << endl;
}