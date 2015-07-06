/*==========================================================
 *
 * MutateAlign.cpp
 *
 * Class for taking a bunch of strands of data, a sequence,
 * and an optional list of mutations, and optimizes the sequence
 * by testing a whole buncha mutations iteratively.
 *
 *========================================================*/

#include "MutateAlign.h"


// define comparison operator heare
bool operator< (const MutScore& ms1, const MutScore& ms2)
{ return ms1.score > ms2.score; }

vector<MutScore> MutateAlign::scoreMutations(const vector<MutInfo>& muts)
{
    // initialize scored output vector
    vector<MutScore> mutscores = vector<MutScore>(muts.begin(),muts.end());
    
    cout << "Scoring";
    cout.flush();
    
    for (int j=0; j<NumEvents; j++)
    {
        // now fill/refill only this alignment, and add the scores from it
        alignments[j].update(this->sequence);
        
        for (int i=0; i<mutscores.size(); i++)
        {
            // sanity check to make sure mutation isn't crappy
            if (muts[i].start > this->sequence.length+4)
                continue;
            // create mutated sequence
            Sequence mutseq(this->sequence,muts[i]);
            // and start the mutated alignments
            mutscores[i].score += alignments[j].scoreMutation(muts[i],mutseq);
        }
        // and now clear this alignment, to save memory
        alignments[j].clear();
        cout << ".";
        cout.flush();
    }
    
    cout << endl;
    cout.flush();
    
    return mutscores;
}

vector<MutScore> MutateAlign::scorePointMutations()
{
    // run through and test all insertions/mutations/deletions and score them
    
    // symbols for insertions/mutations, make static const for persistence
    static uint16_t bases[4] = {'A','C','G','T'};
    
    vector<MutInfo> muts;
    
    for (int i=0; i<sequence.length; i++)
    {
        MutInfo mut;
        mut.start = i;
        
        // deletion
        mut.norig = 1;
        mut.nmut = 0;
        muts.push_back(mut);
        
        // mutations
        mut.norig = 1;
        mut.nmut = 1;
        for (int j=0; j<4; j++)
        {
            // skip over the non-mutation
            if (this->sequence.bases[i] == j)
                continue;
            mut.mut = bases + j;
            muts.push_back(mut);
        }
        
        // insertions
        mut.norig = 0;
        mut.nmut = 1;
        for (int j=0; j<4; j++)
        {
            mut.mut = bases + j;
            muts.push_back(mut);
        }
    }
    
    cout << "Point ";
    
    return scoreMutations(muts);
}

int MutateAlign::testMutations(vector<MutScore> muts)
{
    // how far from a made mutation do we assume scores changed etc?
    const int mutspc = 10;
    
	// how many bases did we mutate?
	int mutbases = 0;
    
    // and sort them by score, descending
    sort(muts.begin(), muts.end());
    // remove all bad (negative) mutations
    while (muts.size() > 0 && muts.back().score < 0)
        muts.pop_back();
    // done if no mutations left
    if (muts.size() == 0)
        return 0;

    cout << "SUPER FAST Testing " << muts.size() << " mutations..." << endl;
    cout.flush();
    
    // vector in which to put ones we didn't use
    vector<MutInfo> mutextra;
    
    // and run through them in order, making them in the sequence
    for (int i=0; i<muts.size(); i++)
    {
        // was it previously invalidated?
        if (muts[i].score < 0)
        {
            mutextra.push_back(MutInfo(muts[i]));
            continue;
        }
        // otherwise, make the mutation
        this->sequence = Sequence(this->sequence,muts[i]);
        // and adjust all future mutations
        cout << "Kept mutation " << i << " at " << muts[i].start 
                << " of " << muts[i].norig << " to " << muts[i].nmut 
                << " with score " << muts[i].score << endl;
        cout.flush();
		// add how many bases we mutated
		mutbases += max(muts[i].norig,muts[i].nmut);

        // we kept a mutation, so we need to modify all subsequent 
        // mutations to account for change
        for (int j=i+1; j<muts.size(); j++)
        {
            // does it overlap with padded, then run it through on the next round instead
            int minind = max(muts[i].start,muts[j].start);
            int maxind = min(muts[i].start+muts[i].nmut,muts[j].start+muts[j].nmut);
            if (minind < maxind+mutspc && muts[j].score > 0)
            {
                // overlaps, do it later
                muts[j].score = -1;
                continue;
            }

            // doesn't overlap, if it starts after, offset the start
            if (muts[j].start >= muts[i].start+muts[i].norig)
                muts[j].start += (muts[i].nmut - muts[i].norig);
        }
    }
    
    if (mutextra.size() > 10)
        mutbases += testMutations(scoreMutations(mutextra));
	
	return mutbases;
}


MutateAlign::MutateAlign(const mxArray* seqarr, const mxArray* evarr,
        AlignParams& par) : params(par), sequence(seqarr)
{
    // how many sequences?
    this->NumEvents = mxGetNumberOfElements(evarr);
    
    // initialize all the events and models
    for (int i=0; i<NumEvents; i++)
        events.push_back(EventData(evarr,i));
    
    // and create the alignments and stuff
    for (int i=0; i<NumEvents; i++)
        alignments.push_back(Alignment(events[i],sequence,params));
}