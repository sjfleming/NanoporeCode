/*==========================================================
 *
 * MutateAlign.h
 *
 * Class for taking a bunch of strands of data, a sequence,
 * and an optional list of mutations, and optimizes the sequence
 * by testing a whole buncha mutations iteratively.
 *
 *========================================================*/

#ifndef _MUTATEALIGN_H_
#define _MUTATEALIGN_H_

#include "mex.h"
#include "EventData.h"
#include "Sequence.h"
#include "Alignment.h"
#include "AlignUtil.h"

class MutateAlign
{
    vector<EventData> events;
    Sequence sequence;
    AlignParams params;
    
    int NumEvents;
    
    vector<Alignment> alignments;
    
public:
    
    MutateAlign(const mxArray* seqarr, const mxArray* evarr, AlignParams& par);
    
    vector<MutScore> scoreMutations(const vector<MutInfo>& muts);
    vector<MutScore> scorePointMutations();
    int testMutations(vector<MutScore> muts);
    
    mxArray* getScores()
    {
        mxArray* arr = mxCreateDoubleMatrix(alignments.size(),1,mxREAL);
        double* pr = mxGetPr(arr);
        
        for (int i=0; i<alignments.size(); i++)
            pr[i] = alignments[i].getMax();
        
        return arr;
    }
    
    mxArray* getSequence()
    {
        mwSize dims[2];
        dims[0] = 1;
        dims[1] = this->sequence.bases.size();
        mxArray* arr = mxCreateCharArray(2,dims);
        uint16_t* pr = (uint16_t*)mxGetPr(arr);
        
        char *str = "ACGT";
        
        for (int i=0; i<dims[1]; i++)
            pr[i] = str[this->sequence.bases[i]];
        
        return arr;
    }
    
    mxArray* getEvents(const mxArray* evarr)
    {
		return EventData::toArray(events, evarr);
    }
};

#endif