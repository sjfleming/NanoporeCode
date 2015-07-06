/*==========================================================
 *
 * viterbi2d.cpp
 *
 * Does simultaneous Viterbi basecalling and alignment on a strand and
 * its complement (or two distinct strands)
 *
 * Arguments: seq1, seq2, model1, model2
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

const int N_BLOCKS = 80000;


// skip and stay probabilities, for both strands
const double skip_prob[2] = {0.13, 0.13};
const double stay_prob[2] = {0.02, 0.02}; // probability of being in a 'stay' state
const double stay_cont[2] = {0.30, 0.30}; // probability of remaining in said state
// maximum number of stays to check for, backtracking 
const int max_stays = 6;

// DIRECTIONALITY NOTES:
// i0 and data0 direction is rows/down
// i1 and data1 direction is cols/right
// thus, L means i1-1, U means i0-1
// so for a step to the right, we look at d1's new level


// a single likelihood element, with overloaded operators
// we need to support addition and max-assignment
struct V2_LIK
{
    double lik;
    double weight;
    int16_t backptr;
    
    V2_LIK() : lik(-inf), weight(1), backptr(-1)
    {}
    V2_LIK(double l) : weight(1), backptr(-1)
    {}
    V2_LIK(double l, double w) : lik(l), weight(w), backptr(-1)
    {}
    V2_LIK(double l, double w, int16_t b) : lik(l), weight(w), backptr(b)
    {}
    
    inline double avg() const
    { return lik/weight; }
};

struct V2_LIKACCUM : public V2_LIK
{
    double avg;
    
    // have only constructor be a copy constructor
    V2_LIKACCUM(const V2_LIK& lik) : V2_LIK(lik)
    {
        avg = lik.avg();
    }
    
    // or one with avg given
    V2_LIKACCUM(const V2_LIK& lik, double av) : V2_LIK(lik)
    {
        avg = av;
    }
};

// add likelihoods and weights
// take backpointer from LHS, but assuming it doesn't really matter
V2_LIK operator+ (const V2_LIK& lik1, const V2_LIK& lik2)
{
    return V2_LIK(lik1.lik+lik2.lik,lik1.weight+lik2.weight,lik1.backptr);
}

V2_LIK& operator+= (V2_LIK& lik1, const V2_LIK& lik2)
{
    lik1 = lik1 + lik2;
    return lik1;
}

// assign max of the two to LHS, comparing averages
V2_LIK& operator<< (V2_LIK& lik1, const V2_LIK& lik2)
{
    if (lik2.avg() > lik1.avg())
        lik1 = lik2;
    return lik1;
}

// or in this case, use V2_LIKACCUM
V2_LIKACCUM& operator<< (V2_LIKACCUM& lik1, const V2_LIK& lik2)
{
    double av2 = lik2.avg();
    if (av2 > lik1.avg)
        lik1 = V2_LIKACCUM(lik2,av2);
    return lik1;
}

// comparison operator
bool operator> (const V2_LIK& lik1, const V2_LIK& lik2)
{ return lik1.avg() > lik2.avg(); }

bool operator< (const V2_LIK& lik1, const V2_LIK& lik2)
{ return lik1.avg() < lik2.avg(); }


// a "block" containing Viterbi model data, back-pointers, and others
class V2_BLOCK
{
public:
    
    V2_LIK state_lik[N_STATES];     // per-state log likelihoods and weights
    V2_BLOCK* neighbors[6];         // pointers to block neighbors, in order {U,L,UL} and their mirrors
    V2_LIK best_lik;                // best log likelihood in block
    int best_state;                 // state of above
    int ind[2];                     // row and col position
    bool prot;                      // is this guy protected/deleteable?
    // do we want some measure of total length to get to this block?
    
    void setLiks(double val)
    {
        for (int i=0; i<N_STATES; i++)
            this->state_lik[i] = V2_LIK(val);
    }
    
    
    V2_BLOCK(int i0, int i1, V2_BLOCK* U, V2_BLOCK* L, V2_BLOCK* UL)
    {
        this->ind[0] = i0;
        this->ind[1] = i1;
        
        // set our private pointers
        for (int i=0; i<6; i++)
            this->neighbors[i] = 0;
        
        this->neighbors[0] = U;
        this->neighbors[1] = L;
        this->neighbors[2] = UL;

        // and DLL pointers on neighbors, if cansies
        if (U) U->neighbors[3] = this;
        if (L) L->neighbors[4] = this;
        if (UL) UL->neighbors[5] = this;
        
        this->best_state = -1;
        this->prot = false;
    }

    ~V2_BLOCK()
    {
        // internal memory is all fine
        // just patch up DLL pointers
        // and DLL pointers on neighbors
        for (int i=0; i<3; i++)
        {
            if (this->neighbors[i]) this->neighbors[i]->neighbors[i+3] = 0;
            if (this->neighbors[i+3]) this->neighbors[i+3]->neighbors[i] = 0;
        }
    }
};

// a class that handles solving stuff
// with internal thingies that do stuff
class V2_SOLVER
{
public:
    // pointers to struct internals
    double* data[2];
    double* datasd[2];
    
    // and how many elements
    int datan[2];
    
    // and model internals
    double* model[2];
    double* modeldev[2];
    double* modelsd[2];
    double* modelsddev[2];
    double* modelwt[2];
    
    // the transition matrices
    double* T[3];
    
    // complement state indices
    int complement[N_STATES];
    
    // and block vector
    vector<V2_BLOCK*> blocks;
    
    
    // populate with a single observation data only, seed row/col
    void populateObs(V2_BLOCK* block, int dataind);
    
    // populate fully from all existing neighbors
    void populateBlock(V2_BLOCK* block);
    
    // get the block avg. likelihoods in a nice output matrix
    mxArray* getBlocks();
    
    // remove worst 10% of blocks
    void trimBlocks();
    
    
    // get a single observation probability, at the ind_th position
    // in dataind dimension
    inline double prob_obs(int dataind, int ind, int state)
    {
        // flip complement strand's observation sequence and states
        if (dataind == 1)
        {
            ind = this->datan[1]-ind-1;
            state = this->complement[state];
        }
        
        return normpdf(this->data[dataind][ind],this->model[dataind][state],this->modeldev[dataind][state])
                *normpdf(this->datasd[dataind][ind],this->modelsd[dataind][state],this->modelsddev[dataind][state]);
    }
    
    
    // constructor, with structs as input directly from mxArray
    V2_SOLVER(const mxArray** args)
    {
        for (int i=0; i<2; i++)
        {
            // get poiners to level vectors
            this->data[i] = mxGetPr(mxGetField(args[i],0,"mean"));
            this->datasd[i] = mxGetPr(mxGetField(args[i],0,"stdv"));
            this->datan[i] = mxGetNumberOfElements(mxGetField(args[i],0,"mean"));
            
            // and model struct fields
            this->model[i] = mxGetPr(mxGetField(args[i+2],0,"level_mean"));
            this->modeldev[i] = mxGetPr(mxGetField(args[i+2],0,"level_stdv"));
            this->modelsd[i] = mxGetPr(mxGetField(args[i+2],0,"sd_mean"));
            this->modelsddev[i] = mxGetPr(mxGetField(args[i+2],0,"sd_stdv"));
            this->modelwt[i] = mxGetPr(mxGetField(args[i+2],0,"weight"));
        }
        
        mexPrintf("Got arrays of size %d and %d\n",this->datan[0],this->datan[1]);
        
        // create a transition matrix
        for (int i=0; i<3; i++)
            this->T[i] = new double[N_STATES*N_STATES];
        
        buildT(T[0], skip_prob[0], stay_prob[0]/(1-stay_cont[0]));
        buildT(T[1], skip_prob[1], stay_prob[1]/(1-stay_cont[1]));
        buildT(T[2], skip_prob[0]*skip_prob[1], stay_prob[0]/(1-stay_cont[0])*stay_prob[1]/(1-stay_cont[1]));
        
        // now take the log of all the transition matrices
        for (int i=0; i<3; i++)
            for (int j=0; j<N_STATES*N_STATES; j++)
                T[i][j] = log(1/inf + T[i][j]);
        
        mexPrintf("Transition matrices built.\n");
        
        // populate complement states
        for (int i=0; i<N_STATES; i++)
            this->complement[i] = complement_state(i);
    }
    
    ~V2_SOLVER()
    {
        // clean up all dem blocks and transition matrix
        for (int i=0; i<this->blocks.size(); i++)
            delete this->blocks[i];
    
        for (int i=0; i<3; i++)
            delete[] this->T[i];
    }
    
    // find a block in blocks with given indices
    V2_BLOCK* findBlock(int i0, int i1)
    {
        for (int i=0; i<this->blocks.size(); i++)
        {
            if (this->blocks[i]->ind[0] == i0 && this->blocks[i]->ind[1] == i1)
                return this->blocks[i];
        }
        return 0;
    }

    
    V2_BLOCK* addBlock(int i0, int i1)
    {
        // first, check if we need to prune blocks
        if (this->blocks.size() >= N_BLOCKS)
            this->trimBlocks();
        
        V2_BLOCK* U = this->findBlock(i0-1,i1);
        V2_BLOCK* L = this->findBlock(i0,i1-1);
        V2_BLOCK* UL = this->findBlock(i0-1,i1-1);

        if (U==0 || UL==0 || L==0)
            mexPrintf("Block not found!\n");
        
        // no blocks found, not origin, don't add
        if (U==0 && UL == 0 && L == 0 && i0 > 0 && i1 > 0)
            return 0;
                
        V2_BLOCK* block = new V2_BLOCK(i0,i1,U,L,UL);
        this->blocks.push_back(block);
        
        // and auto-populate it, based on index
        if (i0 == -1 && i1 == -1)
        {
            // origin block, initialize to zero and set protected
            block->setLiks(0);
            block->prot = true;
            mexPrintf("Origin block.\n");
        }
        else if (i0 == -1)
        {
            // along the top of the matrix
            this->populateObs(block,1);
            block->prot = true;
        }
        else if (i1 == -1)
        {
            // along the left of the matrix
            this->populateObs(block,0);
            block->prot = true;
        }
        else
        {
            // normal block, populate normally
            this->populateBlock(block);
        }
        
        return block;
    }
};

mxArray* backtrace(V2_BLOCK* block, bool setprot);

// sorting operator, for sorting the vector list
// return if a has a higher likelihood, or if it is protected
// so worst blocks are at the end of the list, and we can pop those
bool compBlocks(const V2_BLOCK* a, const V2_BLOCK* b)
{ return (a->best_lik > b->best_lik) || a->prot; }


// fills a block with viterbi likelihoods
void V2_SOLVER::populateBlock(V2_BLOCK* block)
{
    double l_obs[3][N_STATES];
    
    // get observation probabilities for both states
    // even though we may not end up using them....
    
    // calculate as probs first
    // with i=1 being product of both
    for (int j=0; j<N_STATES; j++)
        l_obs[0][j] = log(1/inf + this->prob_obs(0,block->ind[0],j));
    for (int j=0; j<N_STATES; j++)
        l_obs[1][j] = log(1/inf + this->prob_obs(1,block->ind[1],j));
    for (int j=0; j<N_STATES; j++)
    {
        double w1 = this->modelwt[0][this->complement[j]];
        double w2 = this->modelwt[1][j];
        l_obs[2][j] = 2*(l_obs[0][j]*w1+l_obs[1][j]*w2)/(w1+w2);;
    }

    // i==0 is a step from above
    // i==1 is a step from the left
    // i==2 is a step from diagonal
    // (same order as neighbors, *and* arrays!)
    
    for (int i=0; i<3; i++)
    {
        V2_BLOCK* a = block->neighbors[i];
        
        // we don't have a neighbor here
        if (a == 0)
            continue;
        
        // loop through transition matrix, implicitly
        // once for each destination state
        for (int curst=0; curst<N_STATES; curst++)
        {
            double* Tcol = this->T[i]+curst*N_STATES;
            
            // find the maximum likelihood origin state, and keep only that
            // (keep this as local instead of ref to main memory)
            V2_LIKACCUM maxv(block->state_lik[curst]);
            
            double lbase = 0.0;
            if (i<2) lbase = log(skip_prob[i]);
            double curwt = (i<2)?1:2;
            
            for (int oldst=0; oldst<N_STATES; oldst++)
            {
                // speed up & skip
                if (Tcol[oldst] < -100)
                    continue;
                
                V2_LIK lik(Tcol[oldst] + l_obs[i][curst] + lbase,
                        curwt,oldst+(i<<10));
                lik += a->state_lik[oldst];
                // and save it to maxv
                maxv << lik;
            }
            
            // now check for stays, add relevant stay probs
            // (which is both, for diagonal steps)
            lbase = 0.0;
            if (i!=1) lbase += log(stay_prob[0]);
            if (i>0)  lbase += log(stay_prob[1]);
            
            // add current observation probability
            lbase += l_obs[i][curst];
            
            V2_BLOCK* b = a;
            
            for (int nback=0; nback<max_stays; nback++)
            {
                // add accumulated likelihood change (lbase)
                // to backtracked block's likelihood
                V2_LIK lik(lbase, (i<2)?(nback+1):2*(nback+1), curst + (i<<10) + (nback<<12));
                // and save to maxv
                lik += b->state_lik[curst];
                maxv << lik;
                // after each step back, add stay continuation probabilities
                if (i!=1) lbase += log(stay_cont[0]);
                if (i>0)  lbase += log(stay_cont[1]);
                // and observation probabilities
                if (i<2)
                {
                    lbase += log(1/inf + prob_obs(i,b->ind[i],curst));
                }
                else
                {
                    // diagonal, add weighted averages
                    double l1 = log(1/inf + prob_obs(0,b->ind[0],curst));
                    double l2 = log(1/inf + prob_obs(1,b->ind[1],curst));
                    double w1 = this->modelwt[0][this->complement[curst]];
                    double w2 = this->modelwt[1][curst];
                    lbase += 2*(l1*w1+l2*w2)/(w1+w2);
                }
                // and step the block
                b = b->neighbors[i];
                if (b == 0)
                    break;
            }
            
            // now save back to block likelihood
            block->state_lik[curst] = maxv;

            // do we update our block's maximum likelihood value?
            if (maxv > block->best_lik)
            {
                block->best_lik = maxv;
                block->best_state = curst;
            }
        }
    }
    
    if (block->ind[1] == 0)
    {
        mexPrintf("Populated (%d, %d): %0.2f\n",block->ind[0],block->ind[1],block->best_lik.avg());
        cout.flush();
    }
}


void V2_SOLVER::populateObs(V2_BLOCK* block, int dataind)
{
    double p_obs[N_STATES];

    // calculate probabilities
    for (int j=0; j<N_STATES; j++)
        p_obs[j] = this->prob_obs(dataind,block->ind[dataind],j);
   
    // and log-lik them
    for (int j=0; j<N_STATES; j++)
    {
        double l = log(1/inf + p_obs[j]);
        // include a falloff away from the start
        block->state_lik[j] = V2_LIK(l - block->ind[dataind], 1);
        
        if (l > block->best_lik)
        {
            block->best_lik = l;
            block->best_state = j;
        }
    }
    
    mexPrintf("Obs block: %f\n",block->best_lik.avg());
}



// get rid of worst 10% of blocks, or thereabouts
void V2_SOLVER::trimBlocks()
{
    // first, check if we need to do anything
    if (this->blocks.size() < N_BLOCKS)
        return;
    return;
    // run a backtrace to flag last few blocks as protected
    // (from the last 30 added)
    for (int i=1; i<30; i++)
    {
        V2_BLOCK* block = *(this->blocks.end()-i);
        backtrace(block,true);
    }
    
    
    // then sort them all
    sort(this->blocks.begin(), this->blocks.end(), &compBlocks);
    // and get rid of the ones in the back
    for (int i=0; i<N_BLOCKS/10; i++)
    {
        delete this->blocks.back();
        blocks.pop_back();
    }
}



// backtrace optimal path from specified starting point,
// return corresponding mxArray, or just flag as protected
mxArray* backtrace(V2_BLOCK* block, bool setprot)
{
    
    vector<double> liks; // step-by-step likelihood
    vector<int> states;  // state at this step
    vector<int> i0;      // indices of block
    vector<int> i1;      // indices of block
    
    if (block == 0)
        return mxCreateDoubleMatrix(0,0,mxREAL);
    
    int curstate = block->best_state;

    while (block)
    {
        // save values
        liks.push_back(block->state_lik[curstate].avg());
        // 1-based states for matlab
        states.push_back(curstate + 1);
        i0.push_back(block->ind[0]);
        i1.push_back(block->ind[1]);
        // get where we came from       
        int backp = block->state_lik[curstate].backptr;
        // are we done?
        if (backp < 0)
            break;
        // which neighbor direction did we come from?
        int nbor = (backp >> 10) & 3;
        // and how far/how many blocks in that direction?
        // (nback = 0 is adjacent block)
        int nback = (backp >> 12);
        // backtrack the state
        curstate = (backp & (N_STATES-1));
        // and the block
        for (int i=0; i<nback+1; i++)
        {
            if (setprot)
                block->prot = true;
            
            block = block->neighbors[nbor];
        }
    }
    
    if (setprot)
        return 0;
    
    int n = liks.size();
    
    // create output matrix
    mxArray* dpath = mxCreateDoubleMatrix(n,4,mxREAL);
    // and put into matrix, backwards, by columns
    double* arr = mxGetPr(dpath);
    
    for (int i=0; i<n; i++)
    {
        int ri = n-i-1;
        arr[i] = i0[ri];
        arr[i+n] = i1[ri];
        arr[i+2*n] = states[ri];
        arr[i+3*n] = liks[ri];
    }
    
    return dpath;
}

mxArray* V2_SOLVER::getBlocks()
{
    int maxi[2] = {0,0};
    // find out max block indices
    for (int i=0; i<this->blocks.size(); i++)
    {
        for (int j=0; j<2; j++)
            maxi[j] = max(this->blocks[i]->ind[j],maxi[j]);
    }
    
    // add 1 to get dims
    maxi[0]++;
    maxi[1]++;
    
    mxArray* arr = mxCreateDoubleMatrix(maxi[0],maxi[1],mxREAL);
    double *pr = mxGetPr(arr);
    
    // first fill array with all nans
    double nan = mxGetNaN();
    for (int i=0; i<maxi[0]*maxi[1]; i++)
        pr[i] = nan;
    
    for (int i=0; i<this->blocks.size(); i++)
    {
        if (this->blocks[i]->ind[0] < 0 || this->blocks[i]->ind[1] < 0)
            continue;
        
        int ind = maxi[0]*this->blocks[i]->ind[1] + this->blocks[i]->ind[0];
        pr[ind] = this->blocks[i]->best_lik.avg();
    }
    
    return arr;
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // create the solver class
    V2_SOLVER v2s(prhs);
    
    // create top-left origin block, start with zero likelihood
    v2s.addBlock(-1,-1);

    // how many to densely fill, to start
    const int nstart = 75;
    // fill in the edges there
    for (int i=0; i<nstart; i++)
    {
        v2s.addBlock(-1,i);
        v2s.addBlock(i,-1);
    }
    
    mexPrintf("Initialized, starting dense iterations.\n");
    
    // seed the algorithm by filling in some diagonals
    for (int i=0; i<nstart; i++)
        for (int j=0; j<=i; j++)
            v2s.addBlock(j,i-j);
    
    // find best diagonal block
    V2_BLOCK* diag_best = v2s.findBlock(nstart-1,0);
    for (int i=0; i<nstart-1; i++)
    {
        V2_BLOCK* b = v2s.findBlock(i,nstart-1-i);
        if (b->best_lik > diag_best->best_lik)
            diag_best = b;
    }
    
    srand(57);
    
    mexPrintf("Done, starting sparse iterations.\n");
    int bestind = diag_best->ind[0];
    
    for (int idiag=nstart; idiag<2000; idiag++)
    {
        int newind = diag_best->ind[0];
        
        // calculate index of the exact diagonal
        double midind = idiag;
        midind = midind / (1.0 + (double)v2s.datan[1]/(double)v2s.datan[0]);
        // probability that we override the pulled direction, based on current center?
        double p = (bestind-midind)/30.0;
        p = exp(-0.5*p*p);
        
        double d = (double)rand() / RAND_MAX;
        
        // step to the right 
        // a) if we're close to the middle and the best block is to the right
        // b) if we're far from the middle and we are to the left of it
        if ((d < p && newind > bestind) || (d > p && midind > bestind))
            bestind++;
        
        mexPrintf("Sparse iter %d, total size %d, offset %d\n",idiag,v2s.blocks.size(),(int)floor(bestind-midind));
        cout.flush();

        
        V2_BLOCK* newbest = 0;
        
        const int delt = 15;
        
        for (int i=bestind-delt; i<=bestind+delt; i++)
        {
            V2_BLOCK* b = v2s.addBlock(i,idiag-i);
            if (b==0) continue;
            if (newbest == 0 || b->best_lik > newbest->best_lik)
                newbest = b;
        }
        diag_best = newbest;
    }
    
    // now output the scores matrix
    plhs[0] = v2s.getBlocks();
    
    if (nlhs < 2)
        return;
    
    plhs[1] = backtrace(diag_best,false);
    
    // save backtraces from all positions
    /*int dims[2];
    dims[0] = nstart;
    dims[1] = nstart;
    
    plhs[1] = mxCreateCellArray(2,dims);
            
    for (int i=1; i<nstart; i++)
    {
        for (int j=1; j<nstart; j++)
        {
            V2_BLOCK* start = v2s.findBlock(i,j);
            mxSetCell(plhs[1],j*nstart+i,backtrace(start));
        }
    }
    
    if (nlhs < 3)
        return;
    */
    // save full block array of likelihoods and back pointers
    /*plhs[2] = mxCreateCellArray(2,dims);
    
    for (int i=0; i<nb; i++)
    {
        for (int j=0; j<nb; j++)
        {
            // save each block as N_STATESx2 of liks, backpointers
            mxArray* arrOut = mxCreateDoubleMatrix(N_STATES,2,mxREAL);
            mxSetCell(plhs[2],j*nb+i,arrOut);
            V2_BLOCK* block = v2s.findBlock(i,j);
            double* pr = mxGetPr(arrOut);
            for (int k=0;k<N_STATES;k++)
            {
                pr[k] = block->state_lik[k];
                pr[k+N_STATES] = block->backptr[k];
            }
        }
    }*/
}
