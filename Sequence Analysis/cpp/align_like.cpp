/*==========================================================
 *
 * align_like.cpp
 *
 * Aligns an observed sequence to a given reference sequence and model,
 * calculating the likelihood along the way
 *
 * Arguments: Event, Reference
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

// how many likelihood matrices there are
const int AL_LIKS=2;
// how many backpointer possibilities we have
const int AL_BACKS=5;
// and their values/order, for switches
enum AL_MOVE
{
    L_SKIP = 0,
    UL_MATCH = 1,
    U_INSERT = 2,
    U_STAY = 3,
    U_EXTEND = 4,
    IMPLICIT = 255
};

// struct for saving a partial column of a sparse matrix only
struct AL_VEC
{
    int i0;         // starting row of the data
    int i1;         // column of the data
    int length;     // how many data points saved?
    double* liks;   // pointer to the scores array, M x NUM
    uint8_t* steps;     // and the backpointers
    
    // create an empty one
    AL_VEC(int num, int ind0, int ind1) : length(num), i0(ind0), i1(ind1)
    {
        this->liks = (double*)mxCalloc(this->length*AL_LIKS,sizeof(double));
        this->steps = (uint8_t*)mxCalloc(this->length*AL_LIKS,sizeof(uint8_t));
    }
    
    // get a pointer to an element in likelihood array
    double* getPointer(int ind, int arr)
    {
        return this->liks + arr*this->length + (ind-this->i0);
    }
    
    // get a pointer to an element in likelihood array
    uint8_t* getStep(int ind, int arr)
    {
        return this->steps + arr*this->length + (ind-this->i0);
    }
    
    ~AL_VEC()
    {
        mxFree(liks);
        mxFree(steps);
    }
};

typedef unique_ptr<AL_VEC> AL_PTR;

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // probabilities and whatnot, defaults
    double lik_stay = log(0.05);
    double lik_extend = log(0.05);
    double lik_skip = log(0.10);
    double lik_insert = log(0.001);
    // how much to add to observation likelihoods, controls
    // the amount of smith-waterman style alignment
    // (smaller is stricter)
    double lik_offset = 4.5;
    
    // stripe width going up-down, full-frame by default
    int delt = 0;
    
    // read in params, if given
    if (nrhs > 3)
    {
        mxArray* arr;
        arr = mxGetField(prhs[3],0,"stay_prob");
        if (arr) lik_stay = log(mxGetScalar(arr));
        arr = mxGetField(prhs[3],0,"extend_prob");
        if (arr) lik_extend = log(mxGetScalar(arr));
        arr = mxGetField(prhs[3],0,"skip_prob");
        if (arr) lik_skip = log(mxGetScalar(arr));
        arr = mxGetField(prhs[3],0,"insert_prob");
        if (arr) lik_insert = log(mxGetScalar(arr));
        arr = mxGetField(prhs[3],0,"stripe_width");
        if (arr) delt = mxGetScalar(arr);
        arr = mxGetField(prhs[3],0,"lik_offset");
        if (arr) lik_offset = mxGetScalar(arr);
    }
    
    // get poiners to level vector
    EventData event = EventData(prhs[0]);
    
    // and size
    int n0 = event.length;
    int n1 = mxGetNumberOfElements(prhs[2]);
        
    // and get the reference sequence states (ints stored as doubles)
    double* refseq = mxGetPr(prhs[2]);
    
    // make the vector of alignment scores/pointers
    vector<AL_PTR> scores;
    // and create a blank row
    scores.push_back(AL_PTR(new AL_VEC(n0+1,0,0)));
    
    // calculate two-point line coeffs
    double align_m = 0;
    double align_b = 0;
        
    // does it have an existing ref_align? then use that along with
    // a smaller stripe width or something
    if (event.refalend > 0)
    {
        align_m = (event.refalendind-event.refalstartind)/(double)(event.refalend-event.refalstart);
        align_b = event.refalstartind - align_m*event.refalstart;
    }
    else
    {
        // override
        delt = 0;
    }
    
    // full-frame, if not specified
    if (delt == 0)
    {
        delt = n0+1;
        //cout << "Doing full frame alignment..." << endl;
        //cout.flush();
    }

    double maxScore = 0;
    int maxI = 0;
    int maxJ = 0;
    
    // allocate likelihood obs array, full-sized
    // (even if we don't populate it all)
    double* lik_obs = new double[n0+1];

    
    /**************** NOTE NOTE NOTE NOTE NOTE ON INDEXING!!!!!
     * i and j (and all variants thereof) refer to the actual index in an
     * implicit matrix of size (n0+1)*(n1+1), where the first row and
     * column are all 0. Why? Because I said so.
     */
    
    // fill up likelihood matrix
    // (in shorter vertical strips)
    int imid = 1;
    
    for (int j=1; j<=n1; j++)
    {
        // subtract 1 for matlab 1-based indexing
        int refstate = (int)refseq[j-1] - 1;
        
        if (j < event.refalstart || j > event.refalend)
        {
            // calculate imid from linear extrapolation
        	imid = (int)floor(align_m*(j-1) + align_b);
        }
        else
        {
            // use existing refals
            vector<int> inds = event.getrefstates(j);
            if (inds.size() > 0)
                imid = inds.back();
        }
        
        // check all limits and whatnot
        if (imid < 1) imid = 1;
        if (imid > n0) imid = n0;
        
        // upper and lower index to save, inclusive
        int i0 = imid-delt;
        int i1 = imid+delt;
        
        if (i0 < 1) i0 = 1;
        if (i1 > n0) i1 = n0;
        
        // now we know how much we're gonna fill in, allocate corresponding
        // scores vector and pull out their pointers
        scores.push_back(AL_PTR(new AL_VEC(i1-i0+1,i0,j)));
        AL_PTR& curscore = scores[j];
        AL_PTR& prevscore = scores[j-1];
        
        // write all the relevant observation likelihoods
        for (int i=i0; i<=i1; i++)
        {
            lik_obs[i] = lognormpdf(event.mean[i-1],event.model.lev_mean[refstate],event.model.lev_stdv[refstate],event.model.log_lev[refstate])
                + lognormpdf(event.stdv[i-1],event.model.sd_mean[refstate],event.model.sd_stdv[refstate],event.model.log_sd[refstate]);
            lik_obs[i] += lik_offset;
        }
        
        // this is sketchy, because it points to nonexistant memory
        // but, in my defense, i'm, like, _really_ lazy
        double* curlik = curscore->getPointer(0,0);
        uint8_t* curstep = curscore->getStep(0,0);
        double* prevlik = prevscore->getPointer(0,0);
        // and pointers to gap/extension matrix, called "stay" matrix
        double* curstay = curscore->getPointer(0,1);
        uint8_t* curstaystep = curscore->getStep(0,1);
        double* prevstay = prevscore->getPointer(0,1);
        // initialize top element of stay matrix to -inf, not 0
        // since we really super duper can't be there
        curstay[i0] = -inf;
        
        int p0 = prevscore->i0;
        int p1 = prevscore->i0+prevscore->length-1;
        
        for (int i=i0; i<=i1; i++)
        {
            double liks[AL_BACKS] = {0.0,0.0,0.0,-inf,-inf};
            uint8_t likbp[AL_BACKS] = {0,1,2,3,4};
            double lobs = lik_obs[i];
            
            // do implicit zeros if our indices don't work out
            if (i >= p0 && i <= p1)
            {
                liks[L_SKIP] = prevlik[i] + lik_skip;
            }
            else
            {
                liks[L_SKIP] = lik_skip;
                likbp[L_SKIP] = 255;
            }
            
            if (i > p0 && i <= p1)
            {
                liks[UL_MATCH] = prevlik[i-1] + lobs;
            }
            else
            {
                liks[UL_MATCH] = lobs;
                likbp[UL_MATCH] = 255;
            }
            
            // now here is where our stay matrices come into play
            // (no implicit transitions from above; can't stay or extend)
            if (i > i0)
            {
                // start a new stay
                liks[U_STAY] = curlik[i-1] + lobs + lik_stay;
                // or a new stay via insertion
                liks[U_INSERT] = curlik[i-1] + lik_insert;
                // or extend a stay, from stay matrix to stay matrix
                liks[U_EXTEND] = curstay[i-1] + lobs + lik_extend;
            }
            
            // first, update stay matrix with stays or extends
            for (int k=3; k<5; k++)
            {
                if (liks[k] > curstay[i])
                {
                    curstay[i] = liks[k];
                    curstaystep[i] = k;
                }
            }
            // note: stay matrix doesn't go into maxScore
            // just means we don't end on a stay - oh well...

            // ok, now first check transitions from main matrix
            for (int k=0; k<3; k++)
            {
                if (liks[k] > curlik[i])
                {
                    curlik[i] = liks[k];
                    // gotta do this to account for implicit zeros
                    curstep[i] = likbp[k];
                }
            }
            
            // now check "transition" back from stay matrix
            if (curstay[i] > curlik[i])
            {
                curlik[i] = curstay[i];
                curstep[i] = U_STAY;
            }
            
            // and save, if it's the best
            if (curlik[i] > maxScore)
            {
                maxScore = curlik[i];
                maxI = i;
                maxJ = j;
            }
        }
    }
    
    free(lik_obs);
    
    //cout << "Iterations done, starting backtrace..." << endl;
    //cout.flush();
    
    // run backtrace
    // note that this vector will have one value for each event state
    // which is the (one-based!) reference sequence index it aligns with
    // (for stays, we just write the same number multiple times)
    vector<int> inds_i;
    vector<int> inds_j;
    vector<double> ref_like;
    
    int i = maxI;
    int j = maxJ;
    int arr = 0; // which array we're on, main or stay
    
    double score = 1.0;
    
    //cout << "Starting alignment...";
    //cout.flush();
    
    while (i>0)
    {
        //cout << "(" << i << "," << j << "); " << scores.size() << endl;
        //cout.flush();
        
        uint8_t st = *(scores[j]->getStep(i,arr));
        score = *(scores[j]->getPointer(i,arr));
        
        if (score <= 0.0)
            break;
        
        //cout << "(" << i << "," << j << ": " << score << "|" << (int)st << ");" << endl;
        //cout.flush();
        
        switch (st)
        {
            case L_SKIP:
                // skip, don't save
                j--;
                break;
            case UL_MATCH:
                // match, save
                inds_i.push_back(i);
                inds_j.push_back(j);
                ref_like.push_back(score);
                i--;
                j--;
                break;
            case U_INSERT:
                // insert, save but -1 on ref
                inds_i.push_back(i);
                inds_j.push_back(-1);
                ref_like.push_back(score);
                i--;
                break;
            case U_STAY:
                // stay, we're jumping between arrays
                // but only write stuff and change i if it's a stay while
                // in array 1; "stay" in array 0 means we jump to array 1
                // without moving
                if (arr == 1)
                {
                    inds_i.push_back(i);
                    inds_j.push_back(j);
                    ref_like.push_back(score);
                    i--;
                }
                arr = 1-arr;
                break;
            case U_EXTEND:
                // extending stay, save
                inds_i.push_back(i);
                inds_j.push_back(j);
                ref_like.push_back(score);
                i--;
                break;
            case IMPLICIT:
                // we came from an implicit, nonexistant 0
                // (which means skip or match)
                // this is pretty rare, so just quit and leave here
                i = 0;
                break;
            default:
                cout << "BALLS BALLS BALLS" << endl;
                cout.flush();
                i = 0;
                break;
        }
    }
    //cout << "done." << endl;
    //cout.flush();
    
    //cout << "Aligned with " << inds_i.size() << " levels, " << maxScore << endl;
    //cout.flush();
    
    // output the best/last score
    plhs[0] = mxCreateDoubleScalar(maxScore);
    
    if (nlhs > 1)
    {
        // save indices into reference, as copy of event
        plhs[1] = mxDuplicateArray(prhs[0]);
        mxArray* refal = mxGetField(plhs[1],0,"ref_align");
        mxArray* reflik = mxGetField(plhs[1],0,"ref_like");
        double* prefal = mxGetPr(refal);
        double* preflik = mxGetPr(reflik);
        memset(prefal,0,n0*sizeof(double));
        memset(preflik,0,n0*sizeof(double));
        for (int i=0; i<inds_i.size(); i++)
        {
            prefal[inds_i[i]-1] = inds_j[i];
            preflik[inds_i[i]-1] = ref_like[i];
        }
    }    
}
