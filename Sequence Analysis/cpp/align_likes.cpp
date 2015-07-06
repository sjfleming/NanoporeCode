/*==========================================================
 *
 * align_likes.cpp
 *
 * Calculates alignment likelihoods and indices for events and a seq.
 * Also optionally 
 * [scores, events, likes] = align_likes(seq, events, [params])
 *
 *========================================================*/

#include "mex.h"
#include "Alignment.h"

#include <vector>
#include <cstring>


using namespace std;

// Takes event, candidate sequence, and optional params struct

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    AlignParams params; // default parameters
	int numevents = mxGetNumberOfElements(prhs[1]);
	vector<EventData> events;

    // read in params, if given
    if (nrhs > 2)
        params = AlignParams(prhs[2]);

    Sequence sequence(prhs[0]);
		
    // initialize all the events and models
    for (int i=0; i<numevents; i++)
        events.push_back(EventData(prhs[1],i));
		
	// create the output scores array
	plhs[0] = mxCreateDoubleMatrix(numevents,1,mxREAL);
	double* spr = mxGetPr(plhs[0]);
    
    // did we request mapped likelihoods?
	double* lpr = NULL;
    if (nlhs > 2)
    {
        plhs[2] = mxCreateDoubleMatrix(sequence.length+4,1,mxREAL);
        lpr = mxGetPr(plhs[2]);
        memset(lpr,0,(sequence.length+4)*sizeof(double));
    }
        
	// now loop through and do each alignment
	// this in-place updates each event and then
	// destroys the alignment to save memory
	for (int i=0; i<numevents; i++)
	{
		// calculate the matrix
		Alignment align(events[i], sequence, params);
        align.fillColumns();
		// and save the backtrace in the event
		align.backtrace();
		// as well as the output score
		spr[i] = align.getMax();
        // and the likelihood, if requested
        if (lpr > 0)
        {
            // save the aligned likelihoods into lpr
            double lastlik = 0;
            int refind = 1;
            // loop through and find ref_aligns
            for (int j=0; j<events[i].length; j++)
            {
                // save these refinds
                if (events[i].ref_align[j] > 0)
                {
                    // if ref_align[j] == refind, it just updates lastlik
                    for (int k=refind; k<events[i].ref_align[j]; k++)
                        lpr[k+1] += lastlik;
                    lastlik = events[i].ref_like[j];
                    refind = events[i].ref_align[j];
                }
            }
            // and add the last ones on
            for (int k=refind; k<sequence.length+3; k++)
                lpr[k+1] += lastlik;
        }
	}
	
	// and now return the new events too
	if (nlhs > 1)
		plhs[1] = EventData::toArray(events, prhs[1]);
}
