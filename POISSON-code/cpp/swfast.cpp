/*==========================================================
 *
 * swfast.cpp
 *
 * Algorithm for computing a striped Smith-Waterman alignment
 * when the rough alignment is known through other means
 *
 * Arguments: seq0, seq1, paired inds, stripe width
 *
 *========================================================*/

#include "matrix.h"
#include "mex.h"

#include <iostream>
#include <vector>
#include <stdint.h>
#include <cmath>

using namespace std;



const int score_match = 5;
const int score_mismatch = -4;
const int score_insert = -8;


// function in charge of the actual alignment
mxArray* swfast(uint16_t* seq0, uint16_t* seq1, int n0, int n1, double al_m, double al_b, int width)
{
    // m and b are given such that i = m*j + b
    // first, calculate intercepts, so we know where to start and stop
    
    int j0 = (int)floor((-width/2 - al_b)/al_m);
    int j1 = (int)floor((n0+width/2 - al_b)/al_m);
    
    if (j0 < 0) j0 = 0;
    if (j0 >= n1) j0 = n1-1;
    if (j1 < 2) j1 = 2;
    if (j1 > n1) j1 = n1;
    
    // next, allocate arrays, with j0...j1 inclusive
    int* scores = new int[(j1-j0+1)*width]();
    uint8_t* steps = new uint8_t[(j1-j0+1)*width]();
    // and the starting i-indices into above arrays
    int* i0s = new int[(j1-j0+1)];
    // which we then calculate, even if out of range
    for (int j=j0; j<=j1; j++)
        i0s[j-j0] = (int)floor(al_m*j+al_b) - width/2;
    
    // what and where is the max score
    int maxScore = 0;
    int maxI = 0;
    int maxJ = 0;
    
    // go through all relevant columns and calculate
    // (skip the first column to make some if statements easier)
    for (int j=j0+1; j<=j1; j++)
    {
        // so j has range minimum 1 maximum n1 here
        
        // range for the current stripe
        int i0 = i0s[j-j0];
        int i1 = i0 + width - 1;
        // and the stripe to the left
        int p0 = i0s[j-j0-1];
        int p1 = p0 + width - 1;
        
        // clip the values of i0 to 1...n0
        if (i0 < 1) i0 = 1;
        if (i0 > n0) i0 = n0;
        if (i1 < 1) i1 = 1;
        if (i1 > n0) i1 = n0;
        
        // now let's go through and write each output elt.
        // first, create a fake pointer to this column's base
        // and previous column's base
        int* curscore = scores + (j-j0)*width - i0s[j-j0];
        int* prevscore = scores + (j-j0-1)*width - i0s[j-j0-1];
        uint8_t* curstep = steps + (j-j0)*width - i0s[j-j0];
        
        // now loop through i0,i1, which should totes work
        for (int i=i0; i<=i1; i++)
        {
            // best score
            int score=0;
            // best step
            uint8_t step=0;
            
            // note, never gonna have implicit insertions, since those
            // are always negative
            
            // check step from left
            if (i >= p0 && i <= p1)
            {
                int s = prevscore[i] + score_insert;
                if (s > score)
                {
                    score = s;
                    step = 1;
                }
            }
            // and above
            if (i > i0)
            {
                int s = curscore[i-1] + score_insert;
                if (s > score)
                {
                    score = s;
                    step = 2;
                }
            }
            // and above-left, match/mismatch
            if (i > p0 && i <= p1)
            {
                int s = prevscore[i-1] + 
                        ((seq0[i-1]==seq1[j-1])?score_match:score_mismatch);
                if (s >= score)
                {
                    score = s;
                    step = 3;
                }
            }
            else
            {
                // implicit match/mismatch, flag it as such
                int s = (seq0[i-1]==seq1[j-1])?score_match:score_mismatch;
                if (s >= score)
                {
                    score = s;
                    step = 255;
                }
            }
            
            curscore[i] = score;
            curstep[i] = step;
            
            if (score > maxScore)
            {
                maxScore = score;
                maxI = i;
                maxJ = j;
            }
        }
    }
    
    // now run backtrace
    int i = maxI;
    int j = maxJ;
    
    vector<int> inds_i;
    vector<int> inds_j;
    
    while (i > 0 && j > 0)
    {
        int curscore = scores[(j-j0)*width - i0s[j-j0] + i];
        uint8_t curstep = steps[(j-j0)*width - i0s[j-j0] + i];
        
        // hit a 0?
        if (curscore <= 0)
            break;
                
        switch (curstep)
        {
            case 1:
                // step from left, j
                inds_i.push_back(0);
                inds_j.push_back(j);
                j--;
                break;
            case 2:
                // above, i
                inds_i.push_back(i);
                inds_j.push_back(0);
                i--;
                break;
            case 3:
                // les deux
                inds_i.push_back(i);
                inds_j.push_back(j);
                i--;
                j--;
                break;
            case 255:
                // or just done, no previous step
                inds_i.push_back(i);
                inds_j.push_back(j);
                i = 0;
                j = 0;
                break;
            case 0:
            default:
                cerr << "urghhhh bad bad bad day" << endl;
        }
    }

    // and clean up
    delete scores;
    delete steps;
    delete i0s;

    // save to mxarray
    int n = inds_i.size();
    mxArray* arr = mxCreateDoubleMatrix(n,2,mxREAL);
    double *pr = mxGetPr(arr);

    for (int i=0; i<n; i++)
    {
        pr[n-i-1] = inds_i[i];
        pr[2*n-i-1] = inds_j[i];
    }
    
    return arr;
}



/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    if (nrhs < 4)
        mexErrMsgIdAndTxt("CrampBio:swfast","Not enough arguments.");

    // read in sequences
    if (!mxIsChar(prhs[0]) || !mxIsChar(prhs[1]))
        mexErrMsgIdAndTxt("CrampBio:swfast","Arrays not of type char.");
    
    int n0 = mxGetNumberOfElements(prhs[0]);
    int n1 = mxGetNumberOfElements(prhs[1]);
    
    uint16_t* seq0 = (uint16_t*)mxGetPr(prhs[0]);
    uint16_t* seq1 = (uint16_t*)mxGetPr(prhs[1]);
    
    // and the stripe width
    int stripe_width = (int)mxGetScalar(prhs[3]);
    
    // and the index array
    double* pr = (double*)mxGetPr(prhs[2]);
    
    // calculate linear stripe coefficients from given points
    // the order is pr = [i0 j0 i1 j1]
    
    double align_m = (pr[2]-pr[0])/(double)(pr[3]-pr[1]);
	double align_b = pr[0] - align_m*pr[1];

    plhs[0] = swfast(seq0,seq1,n0,n1,align_m,align_b,stripe_width);
}



