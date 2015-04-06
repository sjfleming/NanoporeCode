/*==========================================================
 *
 * Sequence.h
 *
 * Class for manipulating/mutating/whatever a sequence.
 * This implementation is slow, will change if necessary.
 *
 *========================================================*/

#ifndef _SEQUENCE_H_
#define _SEQUENCE_H_

#include "mex.h"
#include "AlignUtil.h"
#include <vector>
#include <iostream>

using namespace std;

struct Sequence
{
    // the length, in states
    int                 length;
    vector<uint8_t>     bases;
    vector<int>         states;
    
    // note that this means state is 4 elements shorter than base
    // which is... kind of necessary
    
    Sequence() : length(0)
    {}
    
    // construct it from a uint16(char) mxArray sequence
    Sequence(const mxArray* arr)
    {
        // make sure it's a char or it's empty
        if (!mxIsChar(arr) && !mxIsEmpty(arr))
            mexErrMsgIdAndTxt("CrampBio:Sequence","Array not of type char.");
        
        this->length = mxGetNumberOfElements(arr);
        uint16_t* pr = (uint16_t*)mxGetPr(arr);
        // and fill the bases
        for (int i=0; i<this->length; i++)
            this->bases.push_back(basetoind(pr[i]));
        
        // the length we save is 4 shorter than the sequence length
        this->length -= 4;
        
        if (this->length > 0)
            this->populateStates();
        else
            this->length = 0;
    }
    
    // construct using a specified mutation
    Sequence(const Sequence& original, const MutInfo& mut)
    {
        // first, copy over the original states
        for (int i=0; i<mut.start; i++)
            this->bases.push_back(original.bases[i]);
        // then insert the new ones
        for (int i=0; i<mut.nmut; i++)
            this->bases.push_back(basetoind(mut.mut[i]));
        // and then the remaining ones
        for (int i=mut.start+mut.norig; i<original.bases.size(); i++)
            this->bases.push_back(original.bases[i]);

        this->length = this->bases.size()-4;
        // then recompute the state indices
        this->populateStates();
    }

private:
    // this function takes the array of bases 0...3 and turns it into
    // an array of states 0...1023
    void populateStates()
    {
        // initialize it with the first few bases
        int curstate = 0;
        for (int i=0; i<4; i++)
            curstate = (curstate << 2) + this->bases[i];
        
        // and then add in the rest
        for (int i=4; i<this->length+4; i++)
        {
            // do we have an invalid state influence?
            if (this->bases[i-4] < 4)
            {
                // no, proceed as normal
                curstate = (N_STATES-1)&((curstate << 2) + this->bases[i]);
                this->states.push_back(curstate);
            }
            else
            {
                // yes, there is a '-' that would influence this state
                curstate = 0;
                this->states.push_back(-1);
            }
        }
    }
};

#endif