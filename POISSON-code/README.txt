

================================= POISSON =================================


Plural Observation-aligning Iterative Sequence-Space Optimization for Nanopores

Contact: tamas@seas.harvard.edu



Code usage notes for obtaining results in manuscript:

 - Main code files start with "poisson_"
 - Path to M13f dataset needs to be set in "poisson_init.m" ("C:\Minion\M13f" by default)
 - For training, path to DNA_CS1 dataset needs to be set in "poisson_train.m"
 - A C++ compiler for mex is required. For windows this means parts of the SDK must be installed.
 - Compilation is performed by typing "make all" in the *Matlab interpreter* (not shell)
 - Compatibility of main code has been tested on both Linux and Windows, Mac is unknown
 - Some analysis/figure code will not run on Linux (align_fasta.m, in particular), but is
   not required to reproduce the results




Executing different parts of the code:


 ---- random trials ----

* poisson_init.m takes a count and trial # and returns a random subset of molecules
  (using trial # as a seed)
* for the manuscript, trials 1-20 were used



 ---- de novo sequencing ----

* Code is in poisson_seq.m
* poisson_seq() starts a single trial, or specify poisson_seq(trial #)
* output is placed in Out/Seq_n.mat and contains self-explanatory arrays
* It might take a while using just a single core
* Results for coverages already found in Out/Seq_n.mat are skipped over, to allow
  paused/continued execution, so if nothing happens, first delete Out/Seq_n.mat



 ---- variant calling ----

* Code is in poisson_variant.m
* Called and outputs identically to poisson_seq, should run without arguments



 ---- model training ----

* Skip and stay params are found in CS_params.conf, and used by above functions
* poisson_train.m loads these params (or defaults, if not found) and iteratively
  optimizes them to maximize the accuracy of the DNA_CS1 dataset
* The function takes no arguments, and attempts to use a few cores, so be warned


 ---- error statistics ----

* poisson_errors.m performs analysis of all Out/Seq_n.mat files to produce figures
  and data shown in supplement section on errors



Outside of this, the purpose and behavior of each individual code file is documented in the files themselves,
and there are a few files not used by the above code but that may be of general utility.

