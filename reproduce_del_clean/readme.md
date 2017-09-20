INSTRUCTIONS TO REPRODUCE SELF-PROPELLED PARTICLE SIMULATIONS AND FIGURES

These simulations are more extensive, and require some cluster time and more complicated setup. (We also provide the key data to generate the figures in the data folder - if you want to use that, skip to step 5.)

Here is our approach:

1) Use the matlab compiler, mcc, to compile the simulation code to an executable. The simulation code is in the folder "simulation_code"\
   You will need to do this either on the cluster you use or a system that can create compatible binaries. 

  mcc -m ensemble_del_psi.m

  This will create an executable with the name "ensemble_del_psi" (assuming you're on *nix)\
  The "ensemble_del_psi" code runs the self-propelled particle simulation amultiple times and averages the results into useful statistics.\
  However, we will need to run this code many different times, over many parameter ranges - I've included some scripts to do this (see below). 

2) Make sure that your cluster has a copy of the Matlab Compiler Runtime (MCR) installed. This can be installed in a local directory in your cluster, so you probably don't need to bug your cluster admin. The version of MCR you install needs to match the version of Matlab you will use to compile the code.

3) Set up scripts to sweep over the parameters you are interested in. In the cluster_tools folder, I've included twosweep.pl, the script I used to generate many different cluster job scripts.

   a) Edit twosweep.pl so that the parameters of interest are being varied. I've included examples for the two figures, twosweep_varykappa.pl and twosweep_varyDpsi.pl
      You will also need to ensure that twosweep.pl points toward the correct location for the MCR, and you might have to change the preamble to match your cluster's architecture.
      See the code's comments for more details!
      
   b) Run twosweep.pl to generate job scripts.

     ./twosweep.pl num_1 num_2 job_name ensemble_del_psi

     will generate num_1xnum_2 scripts named job_name1.sh, job_name2.sh, etc...
     these scripts will vary the first parameter over num_1 points over a range given in twosweep.pl, and the second parameter over num_2 points over the range in twosweep.pl

     To reproduce Fig. 4, use twosweep_varykappa.pl and vary over 10 x 10\
     To reproduce Fig. 5, use twosweep_varyDpsi.pl and vary over 2 x 20

   c) Submit the scripts! I suggest first running only one script, and editing it to have only a small number of iterations (change Nits = 2) to make sure the code runs and saves properly

      The code uses rng('shuffle') to seed the random number generator based on time.
      To avoid having the random numbers generated ever coinciding (and a few other small issues), I put a small delay between each job submission.
      You can do this by editing the script "submitall.sh"  

      While running, the code will generate output in the files t1.out, t2.out, etc...

      With the default parameters, a single run takes about 22 hours on our (fairly old) cluster. 

      Once finished, the code saves a .mat file (size usually tens of MB) with a name set by the savename parameter in twosweep.pl

4) Take the simulation tracks and analyze them. Code for analyzing the tracks is in the analyze folder

    The analysis is mostly centered in the script cellvar_analyze, which generates a matlab structure with analyzed data (od).  The command is

      od = cellvar_analyze(filename_start,Nboots,doplot,doCIerr,nw,nh);

    -filename_start is the first part of the files to be analyzed, i.e. filename_start1.mat, filename_start2.mat, ...\
    -Nboots is the number of bootstrap runs used in estimating the errors in the position-position correlation time (default 50)\
    -if doplot is true, after finishing, the simulation will plot position-position correlations and their exponential fits\
    -doCIerr should be false - it is an attempt to propagate some errors through the CI, but turns out to be both slow and negligible\
    -nw and nh are the number of points simulated for each parameter (i.e. the numbers given to twosweep.pl)

    Unfortunately, this is the opposite order to twosweep.pl - you should enter nw = 20, nh = 2 

    ** You will probably get warnings - for the points where D_psi = 0, the relaxation time is essentially infinite, and the fit fails - we then assume the relaxation time is large (set it = inf)\
    ** This is normal!\
    ** The script will take a few minutes to run

*** Note - to just reproduce the figures, we also provide a data file with typical results from steps 1-4 for both figures in the data folder\
*** These data files, analyzed_vary_kappa.mat and analyzed_vary_Dpsi.mat provide the output of cellvar_analyze.m 

5) Use the scripts to generate figures. Within the folder figures, there are scripts

plot_Dpsi_figure.m and\
plot_kappa_figure.m

Using plot_kappa_figure after loading the od structure for the kappa variation data, and plot_Dpsi_figure after loading the od for the Dpsi variations
will (up to re-arranging the labels) reproduce the figures. (Note plot_kappa_figure.m also shows the full theory lines)

Note - these scripts require the data structure to be named od - also the scripts will change variables in the workspace, etc...

This plotting script also requires the "tight_subplot.m" function (Pekka Kumpulainen), which is included


   

   
