/*==========================================================
 *
 * MutateSequence.cpp
 *
 * Wrapper for C++ code to test given mutations, or point 
 * mutations if none specified.
 * [seq, events, num_mutations] = MutateSequence(seq, events, params, [mutations])
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
	int mutbases = 0;
    if (muts)
        mutbases = mutalign.testMutations(mutalign.scoreMutations(MutInfo::fromArray(muts)));
    else
        mutbases = mutalign.testMutations(mutalign.scorePointMutations());
    
    // return the sequence
    plhs[0] = mutalign.getSequence();
    
    // return the events too
	if (nlhs > 1)
		plhs[1] = mutalign.getEvents(prhs[1]);

    // and the number of mutated bases, if requested
	if (nlhs > 2)
		plhs[2] = mxCreateDoubleScalar(mutbases);
    
}
