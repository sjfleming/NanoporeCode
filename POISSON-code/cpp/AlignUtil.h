/*==========================================================
 *
 * AlignUtil.h
 *
 * Helper function definitions for all the various functions.
 *
 *========================================================*/

#ifndef _ALIGNUTIL_H_
#define _ALIGNUTIL_H_

#include "mex.h"
#include <cmath>
#include <stdint.h>
#include <vector>

using namespace std;


const int N_STATES = 1024;
const double inf=1e300;
#ifndef M_PI
const double M_PI=3.1415926;
#endif
const double log2pi = log(2*M_PI);

// profiling, not presently used
double hightimer();
void ht_prof(int index);

inline int basetoind(uint16_t base)
{
    switch (base)
    {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            return 13;
    }
}

// gaussian dist functions
inline double normpdf(double x, double mu, double sigma)
{
    double d = (x-mu)/sigma;
    d = exp(-0.5*d*d)/sigma/sqrt(2*M_PI);
    return d;
}

inline double lognormpdf(double x, double mu, double sigma, double logsigma)
{
    double d = (x-mu)/sigma;
    return -0.5*(d*d + log2pi) - logsigma;
}

// inverse gaussian dist functions
inline double igpdf(double x, double mu, double lambda)
{
    double d = (x-mu)/mu;
    d = exp(-0.5*d*d*lambda/x);
    return d*sqrt(lambda/(2*M_PI*x*x*x));
}

inline double logigpdf(double x, double mu, double lambda)
{
    double d = (x-mu)/mu;
    return 0.5*(log(lambda/(2*M_PI*x*x*x)) - d*d*lambda/x);
}

// Contains parameters for alignment algorithm
struct AlignParams
{
    // these are the ones that are read by Alignment
    double lik_skip;
    double lik_stay;
    double lik_extend;
    double lik_insert;
    
    // these are the ones loaded from Matlab
    // (this could be cleaner)
    double prob_skip_t;
    double prob_skip_c;
    
    double prob_stay_t;
    double prob_stay_c;
    
    double prob_extend_t;
    double prob_extend_c;
    
    double prob_insert_t;
    double prob_insert_c;

    double lik_offset;
    int stripe_width;
	bool do_verbose;
    
    void setTemplate()
    {
        this->lik_skip = log(this->prob_skip_t);
        this->lik_stay = log(this->prob_stay_t);
        this->lik_extend = log(this->prob_extend_t);
        this->lik_insert = log(this->prob_insert_t);
    }
    
    void setComplement()
    {
        this->lik_skip = log(this->prob_skip_c);
        this->lik_stay = log(this->prob_stay_c);
        this->lik_extend = log(this->prob_extend_c);
        this->lik_insert = log(this->prob_insert_c);
    }
    
    void defaults()
    {
        this->lik_skip = log(0.10);
        this->lik_stay = log(0.05);
        this->lik_extend = log(0.05);
        this->lik_insert = log(0.001);
        
        this->lik_offset = 4.5;
        this->stripe_width = 150;
		this->do_verbose = false;
    }
    
    AlignParams()
    {
        this->defaults();
    }
    
    AlignParams(const mxArray* params)
    {
        this->defaults();
        
        mxArray* arr;
        
        arr = mxGetField(params,0,"skip_t");
        if (arr) this->prob_skip_t = mxGetScalar(arr);
        arr = mxGetField(params,0,"skip_c");
        if (arr) this->prob_skip_c = mxGetScalar(arr);
        arr = mxGetField(params,0,"stay_t");
        if (arr) this->prob_stay_t = mxGetScalar(arr);
        arr = mxGetField(params,0,"stay_c");
        if (arr) this->prob_stay_c = mxGetScalar(arr);
        arr = mxGetField(params,0,"extend_t");
        if (arr) this->prob_extend_t = mxGetScalar(arr);
        arr = mxGetField(params,0,"extend_c");
        if (arr) this->prob_extend_c = mxGetScalar(arr);
        arr = mxGetField(params,0,"insert_t");
        if (arr) this->prob_insert_t = mxGetScalar(arr);
        arr = mxGetField(params,0,"insert_c");
        if (arr) this->prob_insert_c = mxGetScalar(arr);

        arr = mxGetField(params,0,"stripe_width");
        if (arr) this->stripe_width = mxGetScalar(arr);
        arr = mxGetField(params,0,"lik_offset");
        if (arr) this->lik_offset = mxGetScalar(arr);
		arr = mxGetField(params,0,"do_verbose");
        if (arr) this->do_verbose = mxGetScalar(arr)>0;
    }
};


// dumbo little struct to hold mutations
struct MutInfo
{
    int start;
    int norig;
    int nmut;
    uint16_t* mut;

    MutInfo() : start(0), norig(0), nmut(0), mut(0)
    {}
    
    // construct from mxArray
    MutInfo(const mxArray* mutations, int ind)
    {
        start = mxGetScalar(mxGetField(mutations,ind,"start"))-1;
        norig = mxGetNumberOfElements(mxGetField(mutations,ind,"original"));
        nmut = mxGetNumberOfElements(mxGetField(mutations,ind,"mutation"));
        mut = (uint16_t*)mxGetPr(mxGetField(mutations,ind,"mutation"));
    }
    
    // construct from array of mutations
    static vector<MutInfo> fromArray(const mxArray* mutations)
    {
        int nmutations = mxGetNumberOfElements(mutations);
        vector<MutInfo> muts;
        
        for (int i=0; i<nmutations; i++)
            muts.push_back(MutInfo(mutations,i));
        
        return muts;
    }
};

// same as above, but with a score
struct MutScore : MutInfo
{
    double score;

    // initialize with a slightly negative score to get rid of null
    // mutations, which are typically ~+-1e-12    
    MutScore(const MutInfo& mut) : MutInfo(mut), score(-1e-6)
    {}
};

#endif