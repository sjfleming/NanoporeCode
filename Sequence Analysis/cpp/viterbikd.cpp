/*==========================================================
 *
 * viterbikd.cpp
 *
 * Runs Viterbi algorithm on aligned sequences
 *
 * Arguments: levels array
 *
 *========================================================*/

#include "matrix.h"
#include "mex.h"
#include <stdint.h>
#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstring>

#include "viterbi.h"
#include "EventData.h"


using namespace std;


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int nkeep = 0;
    double skip_prob = 0.05;
    double stay_prob = 0.01;
    double mut_min = 0.33;
    double mut_max = 0.75;
    
    // parse params
    if (nrhs > 1)
    {
        mxArray* arr;
        arr = mxGetField(prhs[1],0,"mutations");
        if (arr) nkeep = mxGetScalar(arr);
        arr = mxGetField(prhs[1],0,"skip_prob");
        if (arr) skip_prob = mxGetScalar(arr);
        arr = mxGetField(prhs[1],0,"stay_prob");
        if (arr) stay_prob = mxGetScalar(arr);
        arr = mxGetField(prhs[1],0,"mut_min");
        if (arr) mut_min = mxGetScalar(arr);
        arr = mxGetField(prhs[1],0,"mut_max");
        if (arr) mut_max = mxGetScalar(arr);
    }
     
    // how many sequences?
    int n_seq = mxGetNumberOfElements(prhs[0]);
    
    vector<EventData> events;
    
    // contain the dynamically growing array of scores
    vector<VK_PTR> scores;

    // initialize all the events and models
    for (int i=0; i<n_seq; i++)
        events.push_back(EventData(prhs[0],i));

    // so now all seqinds are at their first nonzero location
    // and refind is set to the first aligned reference state
    
    // let's make a blank, dummy initial likelihood
    scores.push_back(VK_PTR(new VK_LIK()));
    VK_PTR& lik0 = scores.back();
        
    for (int i=0; i<N_STATES; i++)
    {
        lik0->liks[i] = 0;
        lik0->backptrs[i] = -1;
        lik0->fwdprobs[i] = 1.0/N_STATES;
    }
	
	cout << "Viterbi";
    
    double* obs = new double[N_STATES*n_seq];
    
    int refind = events[0].refstart;
    for (int i=0; i<n_seq; i++)
        refind = min(refind,events[i].refstart);
    
    // now loop through each level
    while (1)
    {
        // calculate log-normalized observations
        //double obs[N_STATES*n_seq];
        
        memset(obs,0,N_STATES*sizeof(double));
        
        int nlik = 0;
        for (int k=0; k<n_seq; k++)
        {
            vector<int> inds = events[k].getrefstates(refind);
            
            // does this strand have any states that line up with this ref state?
            if (inds.size() == 0)
                continue;
            
            nlik++;
            
            // now, average all these thingies, and stuff
            double lvl = 0;
            double sd = 0;
            for (int j=0; j<inds.size(); j++)
            {
                lvl += events[k].mean[inds[j]];
                sd += events[k].stdv[inds[j]];
            }
            lvl = lvl/inds.size();
            sd = sd/inds.size();

            for (int j=0; j<N_STATES; j++)
                obs[j*n_seq + nlik-1] = log(1/inf + events[k].model.getProb(lvl,sd,j));
        }
        
        // how many align here?
        int nalhere = 0;
        for (int k=0; k<n_seq; k++)
            if (refind >= events[k].refstart && refind <= events[k].refend)
                nalhere++;
        
        // if too few aligned strands, skip over this state (or end, if done)
        if (nlik <= nalhere*0.2)
        {
            // stop, if no more strands
            if (nalhere == 0)
                break;
            
            // or just skip this one
            refind++;
            continue;
        }
        
        if (nlik > 1)
        {
            // sort obs likelihoods, ascendingally
            for (int j=0; j<N_STATES; j++)
                sort(obs+j*n_seq,obs+j*n_seq+nlik);
            // and then calculate means likelihood of the last few elements only
            int nskip = floor(nlik*0.25);
            if (nskip > nlik-2) nskip = 0;
            for (int j=0; j<N_STATES; j++)
            {
                double lik = 0.0;
                for (int k=nskip; k<nlik; k++)
                    lik += obs[j*n_seq + k];
                // now put it in the first column, this is shitty code
                // but again, suck it
                obs[j] = lik / (nlik - nskip);
            }
        }
        else
        {
            // nothing really to do
            for (int j=0; j<N_STATES; j++)
                obs[j] = obs[j*n_seq];
        }

        // save a new likelihood
        double curskipprob = skip_prob;// + 0.3*(nalhere-nlik)/(double)n_seq;
        double curstayprob = stay_prob;// + 0.1*(nalhere-nlik)/(double)n_seq;
        scores.push_back(VK_PTR(new VK_LIK(scores.back(),obs,curskipprob,curstayprob)));        

        // and advance
        refind++;
        
        if (refind%100 == 0)
        {
            mexPrintf(".");
            cout.flush();
        }
    }
    delete[] obs;
    mexPrintf("\n");
    cout.flush();    
    
    // now output the resulty thingy
    int n = scores.size() - 1;
    plhs[0] = mxCreateDoubleMatrix(n,1 + nkeep,mxREAL);
    double *pr = mxGetPr(plhs[0]);
    
    // run backtrace on main path
    // find largest final likelihood
    double* mlik = max_element(scores.back()->liks,scores.back()->liks+N_STATES);
    int curst = mlik - scores.back()->liks;
    
    for (int i=n-1; i>=0; i--)
    {
        // 1-based states
        pr[i] = curst + 1;
        curst = scores[i+1]->backptrs[curst];
    }
    
    // and now for mutations (if needed)
    if (nkeep == 0)
        return;    
    
	
	
    // create a transition matrix
    double* T = new double[N_STATES*N_STATES];
    
    buildT(T, skip_prob, stay_prob);
        
    // now make sure transition matrix knows about stays properly
    for (int i=0; i<N_STATES; i++)
        T[i*(1+N_STATES)] = stay_prob;

    
    for (int k=0; k<nkeep; k++)
    {
        for (int i=n-1; i>=0; i--)
        {
            // 1-based states
            pr[i + (k+1)*n] = curst + 1;
            curst = scores[i+1]->randbp(curst,mut_min+(mut_max-mut_min)*k/(double)nkeep,T);
        }
    }
        
    delete[] T;
}
