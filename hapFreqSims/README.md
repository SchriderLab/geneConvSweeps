# Custom simulation code for tracking haplotype frequencies during (pseudo-)soft sweeps

This directory contains the code for simulating sweeps with gene conversion or recurrent mutation and measuring the frequencies of each individual haplotype that the beneficial mutation is either copied onto via gene conversion or originates on via recurrent mutation.

## Contents

1. simTestRecMut.py: This is the custom Python simulation that was used to examine the frequencies of distinct haplotypes after a sweep. It was initially used to check the results of the SLiM simulations in the parent directory for a subset of parameters (hence the name), and then expanded to include recurrent mutation and track the distinct haplotypes created by gene conversion or recurrent mutation. It also calculates several summaries of haplotypic diversity immediately after fixation (including the number of distinct haplotypes, H2/H1 in soft sweeps, etc.). The command-line arguments are the following (in order): the (constant) population size, the relative gene conversion rate (i.e. the g parameter from the paper), the relative recurrent mutation rate (the m parameter from the paper), the selection coefficient, the dominance coefficient, the population size rescaling factor Q, and the number of replicates to run. You can run the script directly or use the launch script below:

2. `launchAllGcAndRecMutSims.py`: Code for launching the simulations on a high-performance computing cluster. The code is written assuming that the cluster uses the `SLURM` scheduler and that the desired partition name is `general`, so the code may have to be modified to run on your computing resources. It takes no command-line arguments.

3. `combineAllGcAndRecMutSimResults.py`: Aggregates all results into a new directory called `simOutcombined/`

4. `hap_freq_analysis.ipynb`: The notebook that plots the results.
