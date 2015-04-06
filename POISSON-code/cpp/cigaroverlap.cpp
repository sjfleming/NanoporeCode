/*==========================================================
 *
 * cigaroverlap.cpp
 *
 * Get substring from a CIGAR string, sequence, cigar start index,
 * and start, end indices to ref
 * 
 * Returns the subsequence, padded with '-' where applicable
 *
 *========================================================*/

#include "matrix.h"
#include "mex.h"
#include <stdint.h>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    if (nrhs != 5)
        mexErrMsgIdAndTxt("CrampBio:cigaroverlap","Five inputs required.");
    for (int i=0; i<2; i++)
        if (!mxIsChar(prhs[i]))
            mexErrMsgIdAndTxt("CrampBio:cigaroverlap","2 inputs must be type char!");
    
    uint16_t* cigar = (uint16_t*)mxGetPr(prhs[0]);
    uint16_t* seq = (uint16_t*)mxGetPr(prhs[1]);
    
    int ncigar = mxGetNumberOfElements(prhs[0]);
    int ns = mxGetNumberOfElements(prhs[1]);
    
    int indc = 0; // index in cigar string
    // index in reference, start at given value
    // (less one for indexing reasons)
    int indr = (int)mxGetScalar(prhs[2]) - 1;
    // and the desired start and end ranges
    int refstart = (int)mxGetScalar(prhs[3]) - 1;
    int refend = (int)mxGetScalar(prhs[4]) - 1;
    int inds = 0; // index in sequence
    
    // this is the current integer value of our cigar string
    int cigval = 0;
    
    // and the growing sequence
    vector<uint16_t> subseq;
    
    // pad with leading dashes, if indr > refstart
    while (indr > refstart)
    {
        subseq.push_back('-');
        refstart++;
    }
        
    
    // now loop through cigar string to end (or end of seq/cigar)
    while (indr < refend && inds < ns && indc < ncigar)
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
                    // match or mismatch, step both (and save)
                    for (int i=0; i<cigval; i++)
                    {
                        if (indr >= refstart)
                            subseq.push_back(seq[inds]);
                        indr++;
                        inds++;
                    }
                    break;
                case 'I':
                    // insertion, relative to ref, so advance sequence
                    for (int i=0; i<cigval; i++)
                    {
                        if (indr >= refstart)
                            subseq.push_back(seq[inds]);
                        inds++;
                    }
                    break;
                case 'D':
                    // deletion, relative to ref, so advance ref
                    indr += cigval;
                    break;
                case 'S':
                    cerr << "S operation not implemented!" << endl;
                    break;
                case 'H':
                    // hard clip, do nothing, since we're using the clipped seq
                    break;
            }
            cigval = 0;
        }
        indc++;
    }
    
    // now make sure we have enough bases
    while (indr < refend)
    {
        subseq.push_back('-');
        indr++;
    }
    
    // and save the sequence
    mwSize dims[2];
    dims[0] = 1;
    dims[1] = subseq.size();
    mxArray* arr = mxCreateCharArray(2,dims);
    uint16_t* pr = (uint16_t*)mxGetPr(arr);

    for (int i=0; i<dims[1]; i++)
        pr[i] = subseq[i];

    plhs[0] = arr;
    
}