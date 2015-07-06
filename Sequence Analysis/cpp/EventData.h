/*==========================================================
 *
 * EventData.h
 *
 * Classes to hold event and model data and stuff.
 *
 *========================================================*/

#ifndef _EVENTDATA_H_
#define _EVENTDATA_H_

#include "mex.h"
#include <vector>
#include <algorithm>

#include "AlignUtil.h"
#include "Sequence.h"

using namespace std;

// internal model struct
struct ModelData
{
    // supplied values
    double lev_mean[N_STATES];
    double lev_stdv[N_STATES];
    double sd_mean[N_STATES];
    double sd_stdv[N_STATES];
    
    // computed values
    double log_lev[N_STATES];
    double log_sd[N_STATES];
    // inverse gaussian parameters
    double sd_lambda[N_STATES];
    double log_lambda[N_STATES];
    
    double skip_prob;
    double stay_prob;
    double extend_prob;
    
     // initialize the model from a Matlab array
    void init(const mxArray* arr, int ind)
    {
        double* levmean = mxGetPr(mxGetField(arr,ind,"level_mean"));
        double* levstdv = mxGetPr(mxGetField(arr,ind,"level_stdv"));
        double* sdmean = mxGetPr(mxGetField(arr,ind,"sd_mean"));
        double* sdstdv = mxGetPr(mxGetField(arr,ind,"sd_stdv"));
        
        this->skip_prob = mxGetScalar(mxGetField(arr,ind,"skip_prob"));
        this->stay_prob = mxGetScalar(mxGetField(arr,ind,"stay_prob"));
        this->extend_prob = mxGetScalar(mxGetField(arr,ind,"extend_prob"));

        // copy over and also precalculate logarithmed versions
        for (int i=0; i<N_STATES; i++)
        {
            this->lev_mean[i] = levmean[i];
            this->lev_stdv[i] = levstdv[i];
            this->sd_mean[i] = sdmean[i];
            this->sd_stdv[i] = sdstdv[i];
            this->log_lev[i] = log(levstdv[i]);
            this->log_sd[i] = log(sdstdv[i]);
            
            this->sd_lambda[i] = pow(this->sd_mean[i],3)/pow(this->sd_stdv[i],2);
            this->log_lambda[i] = log(this->sd_lambda[i]);
        }
    }
    
    ModelData(const mxArray* arr, int ind)
    {
        this->init(arr, ind);
    }
    
    ModelData(const mxArray* arr)
    {
        this->init(arr, 0);
    }
    
    ModelData()
    {
    }

    // probability helper function. most likely i shouldn't have this.
    inline double getProb(double lev, double sd, int ind)
    {
        return normpdf(lev,this->lev_mean[ind],this->lev_stdv[ind])
            *normpdf(sd,this->sd_mean[ind],this->sd_stdv[ind]);
    }
};


// struct for an event
struct EventData
{
    // internal model structure corresponding to event
    ModelData model;
    // and the sequence data associated with the event
    // (the 2d or 1d original sequence)
    Sequence sequence;
    
    // how many levels?
    int length;
    // beginning and end of reference
    int refstart;
    int refend;

    vector<double> mean;
    vector<double> stdv;
    vector<double> ref_align;
    vector<double> ref_index;
    vector<double> ref_like;
    
    // populate ref_index densely with indices and such
    void updaterefs()
    {
        double al_m, al_b;  // alignment coefficients to extend past start
        int ra0, ra1;       // indices in ref_align
        
        ra0 = -1;
        ra1 = -1;
        refstart = -1;
        refend = -1;
        // ok, find first and last nonnegative index in ref_align
        for (ra0=0; ra0<this->length; ra0++)
            if (ref_align[ra0] > 0)
                break;
        for (ra1=this->length-1; ra1>=0; ra1--)
            if (ref_align[ra1] > 0)
                break;

        // if we don't have any, abort
        if (ra0 < 0 || ra1 < 0)
        {
            ref_index = vector<double>();
            return;
        }
        
        // now set ref start and end
        refstart = ref_align[ra0];
        refend = ref_align[ra1];
        
        // otherwise, initialize array as a copy
        ref_index = ref_align;
        
        // calculate linear alignment and stuff
        al_m = (ref_align[ra1]-ref_align[ra0])/(double)(ra1-ra0);
        al_b = ref_align[ra0] - al_m*ra0;
        // now fill in all of it
        int lastal = -1; // previous ref_align-populated state
        for (int i=0; i<this->length; i++)
        {
            // before first and last aligned states
            if (i < ra0 || i > ra1)
            {
                // just use linear extrapolation
                ref_index[i] = al_m*i + al_b;
            }
            // we have an alignment
            else if (ref_align[i] > 0)
            {
                if (lastal > 0)
                {
                    // loop from lastal + 1
                    double m = (ref_align[i]-ref_align[lastal])/(i-lastal);
                    // and fill in all the ref_indexes
                    for (int j=lastal+1; j<i; j++)
                        ref_index[j] = m*(j-lastal) + ref_align[lastal];
                }
                lastal = i;
                // ref_index[i] already is ref_align[i]
            }
        }
        
        //for (int i=0; i<length; i++)
        //   mexPrintf("Ind: %d, RefAl: %0.1f, RefInd: %0.1f\n",i,ref_align[i],ref_index[i]);
    }
    
    // find just the one single index
    int getrefstate(int refind)
    {
        // if empty, just play dumb
        if (ref_index.size() == 0)
            return 0;
        // find first position greater than or equal to refind
        vector<double>::iterator it = std::lower_bound(ref_index.begin(),ref_index.end(),refind);
        // if refind is less than ref_index[0], this will be fine, since it == ref_index.begin() then
        // and if refind is more than the end of ref_index, also fine, since default behavior is to return
        // ref_index.end() if not found
        return it - ref_index.begin();
    }

    // find all the indices that lie within a certain ref index
    // (only if it exists, yadda yadda)
    vector<int> getrefstates(int refind)
    {
        vector<int> inds;
        
        // find index of refind in ref_index array
        vector<double>::iterator it = std::find(ref_index.begin(),ref_index.end(),refind);
        // if we didn't find it, return nada
        if (it == ref_index.end())
            return inds;
        
        int i = it - ref_index.begin();
        inds.push_back(i);
        for (i++; i<length && ref_align[i] <= refind; i++)
            if (ref_align[i] > 0)
                inds.push_back(i);
        
        return inds;
    }
    
    
    // return event(s) as a Matlab array
    static void toArray(EventData& event, mxArray* arr, int ind)
    {
        // save indices into reference, as copy of event
        mxArray* refal = mxGetField(arr,ind,"ref_align");
        mxArray* reflik = mxGetField(arr,ind,"ref_like");
        double* prefal = mxGetPr(refal);
        double* preflik = mxGetPr(reflik);
        for (int i=0; i<event.length; i++)
        {
            prefal[i] = event.ref_align[i];
            preflik[i] = event.ref_like[i];
        }
   }
    
    // output a single array back to matlab
    static mxArray* toArray(EventData& event, const mxArray* evarr)
    {
        mxArray* arr = mxDuplicateArray(evarr);
        toArray(event,arr,0);
        return arr;
    }
    
    // output a vector of arrays back to matlab
    static mxArray* toArray(vector<EventData>& events, const mxArray* evarr)
    {
        // make output copy of array
        mxArray* arr = mxDuplicateArray(evarr);
        // and write each event to it
        for (int i=0; i<events.size(); i++)
            toArray(events[i],arr,i);
        // and return that array
        return arr;
    }
    
    void init(const mxArray* arr, int ind)
    {
        double* pr;
        // initialize the big uns
        this->model = ModelData(mxGetField(arr,ind,"model"));
        this->sequence = Sequence(mxGetField(arr,ind,"sequence"));
        // and initialize the internal data stores
        this->length = mxGetNumberOfElements(mxGetField(arr,ind,"mean"));
        pr = mxGetPr(mxGetField(arr,ind,"mean"));
        this->mean = vector<double>(pr,pr+this->length);
        pr = mxGetPr(mxGetField(arr,ind,"stdv"));
        this->stdv = vector<double>(pr,pr+this->length);
        pr = mxGetPr(mxGetField(arr,ind,"ref_align"));
        this->ref_align = vector<double>(pr,pr+this->length);
        pr = mxGetPr(mxGetField(arr,ind,"ref_like"));
        this->ref_like = vector<double>(pr,pr+this->length);

        
        // and populate ref_index to start
        this->updaterefs();
    }
    
    EventData(const mxArray* arr, int ind)
    {
        this->init(arr, ind);
    }
    
    EventData(const mxArray* arr)
    {
        this->init(arr, 0);
    }
    
    EventData() : length(0)
    {}
};


#endif