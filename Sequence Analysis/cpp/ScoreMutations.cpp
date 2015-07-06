/*==========================================================
 *
 * ScoreMutations.cpp
 *
 * Calculates the scores of given mutations (or point mutations)
 * without actually making them.
 * [scores] = ScoreMutations(seq,events,params,[mutations])
 *
 *========================================================*/


#include "mex.h"
#include "Alignment.h"
#include "MutateAlign.h"

using namespace std;

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    AlignParams params;

    // read in params, if given
    if (nrhs > 2)
        params = AlignParams(prhs[2]);
    
    const mxArray* muts = 0;
    // and mutations, if given
    if (nrhs > 3)
        muts = prhs[3];
    
    // create mutatealign class, which seeds everything with default
    // alignment for all the strands
    MutateAlign mutalign(prhs[0],prhs[1],params);
    
    // and now run the mutations
    vector<MutScore> mutscores;
    if (muts)
        mutscores = mutalign.scoreMutations(MutInfo::fromArray(muts));
    else
        mutscores = mutalign.scorePointMutations();
    
    // return the scores
    plhs[0] = mxCreateDoubleMatrix(mutscores.size(),1,mxREAL);
    double* lpr = mxGetPr(plhs[0]);
    for (int i=0; i<mutscores.size(); i++)
        lpr[i] = mutscores[i].score;
}
