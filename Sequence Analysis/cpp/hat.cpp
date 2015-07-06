/* HIGH ACCURACY TIMER
 *
 * compile command (needs windows SDK):
 * mex -O hat.c
 *
 * Ivo Houtzager
 * modified by Tamas Szalay
 */

#ifndef __linux__

#include "windows.h"
#include "mex.h"
#include <stdio.h>

double hightimer()
{
    HANDLE hCurrentProcess = GetCurrentProcess();
    DWORD_PTR dwProcessAffinity;
    DWORD_PTR dwSystemAffinity;    
    static LARGE_INTEGER frequency;
    LARGE_INTEGER counter;
    double sec_per_tick, total_ticks;
    
    /* force thread on first cpu */
    //GetProcessAffinityMask(hCurrentProcess,&dwProcessAffinity,&dwSystemAffinity);
    //SetProcessAffinityMask(hCurrentProcess, 1);
    
	/* retrieves the frequency of the high-resolution performance counter */
    QueryPerformanceFrequency(&frequency);
    
    /* retrieves the current value of the high-resolution performance counter */
    QueryPerformanceCounter(&counter);
    
     /* reset thread */
    //SetProcessAffinityMask(hCurrentProcess,dwProcessAffinity);

	/* time in seconds */
    sec_per_tick = (double)1/(double)frequency.QuadPart;
    total_ticks = (double)counter.QuadPart;  
    return sec_per_tick * total_ticks;
}

void ht_prof(int index)
{
    const int n_prof = 16;
    static double ptime[n_prof] = {0};
    static double pstart[n_prof] = {0};
    static int pcount[n_prof] = {0};
    
    if (index > 0)
    {
        // start a timer
        pstart[index-1] = hightimer();
    }
    else if(index < 0)
    {
        // stop the timer and save
        index = -index;
        double dt = hightimer() - pstart[index-1];
        ptime[index-1] += dt;
        pcount[index-1]++;
    }
    else
    {
        // display
        for (int i=0; i<n_prof; i++)
        {
            if (pcount[i] == 0)
                continue;
            
            mexPrintf("Timer %d: %0.3g/%0.3g\n",i+1,ptime[i]/pcount[i],ptime[i]);
        }
        for (int i=0; i<n_prof; i++)
        {
            pcount[i] = 0;
            ptime[i] = 0;
        }
    }
}

#else

double hightimer() {return 0;}
void ht_prof(int index) {}



#endif

