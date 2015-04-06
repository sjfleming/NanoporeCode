/*==========================================================
 *
 * viterbi.h
 *
 * Helper functions/transition matrix for 1d, 2d C++ Viterbi algorithms
 * The Viterbi code is older code, barely used, and pretty junky...
 *
 *========================================================*/

#ifndef _VITERBI_H_
#define _VITERBI_H_

#include <stdint.h>
#include "mex.h"
#include <vector>

#include "AlignUtil.h"

using namespace std;

// previous state macro/inline function
inline int prev_state(int state, int ind)
{ return (state >> 2)+(ind<<8); }
inline int next_state(int state, int ind)
{ return ((state << 2)&(N_STATES-1))+ind; }
// and with multiple advances
inline int prev_state(int state, int ind, int nsteps)
{ return (state >> (2*nsteps))+(ind<<(10-2*nsteps)); }
inline int next_state(int state, int ind, int nsteps)
{ return ((state << (2*nsteps))&(N_STATES-1))+ind; }

// returns complement of a base
int complement_state(int state)
{
    int comp = 0;
    
    // reverse two bits at a time, taking xor of said bits
    for (int i=0; i<5; i++)
    {
        comp <<= 2;
        comp += (state&3)^3;
        state >>= 2;
    }
    return comp;
}

inline void normvec(double* vec)
{
    double tot = 0;
    for (int i=0; i<N_STATES; i++)
        tot += vec[i];
    tot = 1.0/tot;
    for (int i=0; i<N_STATES; i++)
        vec[i] *= tot;
}


void buildT(double* T, double skip_prob, double stay_prob)
{
    // initialize to 0
    memset(T,0,N_STATES*N_STATES*sizeof(double));
        
    // number of skips to consider ('1 skip' is a normal step)
    const int nskip = 4;
    
    // all destination states
    for (int curst=0; curst<N_STATES; curst++)
    {
        // current row (or column, or whatever)
        // small steps in T represent origin states
        double* Tcol = T + N_STATES*curst;
        
        // a loop for each number of skips
        double sp = 0.25;
        for (int j=1; j<=nskip; j++)
        {
            // loop through all possible new states
            for (int k=0; k < 1<<(2*j); k++)
            {
                int prevst = prev_state(curst, k, j);
                Tcol[prevst] += sp;
            }
            // and change prob for next time through
            sp = sp*0.25*skip_prob;
        }
        
        // stay probabilities are done in the calculation code
        // and are not in the T matrix
    }
    
    // now go through and normalize T, in the other dimension
    // (including the stay probabilities in this calculation)
    /*for (int prevst=0; prevst<N_STATES; prevst++)
    {
        double ptot = stay_prob;
        
        for (int curst=0; curst<N_STATES; curst++)
            ptot += T[curst*N_STATES + prevst];
        
        ptot = 1.0/ptot;
        
        for (int curst=0; curst<N_STATES; curst++)
            T[curst*N_STATES + prevst] *= ptot;
    }*/
}

// holds viterbi likelihoods etc etc
struct VK_LIK
{
    double liks[N_STATES];              // accumulated likelihood so far, best-path
    int backptrs[N_STATES];             // viterbi backpointers
    double fwdprobs[N_STATES];          // forward probability values
    
    VK_LIK(unique_ptr<VK_LIK>& prevlik, double* obs, double skip_prob, double stay_prob);
    VK_LIK() {};
    
    // get a random backpointer, for a given state, using a given
    // transition matrix (not log-form)
    int randbp(int curstate, double atten, double* T);
    
};
typedef unique_ptr<VK_LIK> VK_PTR;

VK_LIK::VK_LIK(VK_PTR& prevlik, double* obs, double skip_prob, double stay_prob)
{
    // number of skips to consider ('1 skip' is a normal step)
    const int nskip = 3;
    const double skip_lik = log(skip_prob);
    const double stay_lik = log(stay_prob);
    
    // all destination states
    for (int curst=0; curst<N_STATES; curst++)
    {
        // best numbers for this state
        double maxlik = -inf;
        int maxptr = -1;
        // forward algorithm probability
        double fwdprob = 0.0;


        // a loop for each number of skips
        double sp = 0.25;
        double lsp = log(0.25);
        
        for (int j=1; j<=nskip; j++)
        {
            // loop through all possible states we came from
            for (int k=0; k < 1<<(2*j); k++)
            {
                int prevst = prev_state(curst, k, j);
                
                double l = obs[curst] + lsp;
                l += prevlik->liks[prevst];
                fwdprob += sp * prevlik->fwdprobs[prevst];
                if (l > maxlik)
                {
                    maxlik = l;
                    maxptr = prevst;
                }
            }
            // and change prob for next time through
            sp = sp*0.25*skip_prob;
            lsp = lsp + log(0.25) + skip_lik;
        }
        
        // now finally check stays, add it to self-transition prob
        double l = obs[curst]+stay_lik+prevlik->liks[curst];
        if (l > maxlik)
        {
            maxlik = l;
            maxptr = curst;
        }

        fwdprob += stay_prob * prevlik->fwdprobs[curst];

        // multiply forward probability by the observation probability
        fwdprob *= exp(obs[curst]);

        // and partially sort the array to get best few likelihoods
        // save the best likelihood only into lik array
        this->liks[curst] = maxlik;
        this->backptrs[curst] = maxptr;
        this->fwdprobs[curst] = fwdprob;
    }
    
    // normalize forward probabilities
    normvec(this->fwdprobs);
}

// generate a random backpointer using forward probs (it ... sort of makes sense?)
int VK_LIK::randbp(int curstate, double atten, double* T)
{
    // random value 0...1
    double r = rand() / (double(RAND_MAX)+1);
    // and value to compare to
    double cumsum = 0;
    // now renormalize probs multiplied by transitions
    double probs[N_STATES];
    
    for (int i=0; i<N_STATES; i++)
        probs[i] = T[i + curstate*N_STATES] * pow(this->fwdprobs[i],atten);

    // and normalize it, booyah
    normvec(probs);
    
    // and find the one that is hit by r
    for (int i=0; i<N_STATES; i++)
    {
        cumsum += probs[i];
        if (r < cumsum)
            return i;
    }
    // if we didn't find one yet, just some numerical slop, it should
    // be the last state
    // (this doesn't happen in practice)
    return N_STATES-1;
}

#endif
