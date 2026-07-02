# Code for misc tests of simualtion accoracy (rescaling, and comparison between SLiM and my custom Python simulator)

This directory contains code for running two comparisons: 1) `SLiM` simulations of pseudo-soft sweeps with different values of the rescaling coefficient (Q), 2) comparing the results between SLiM and the custom Python simulation code found in the ../hapFreqSims directory.

## Contents

1. `launchComparisonPySims.py`: Code for launching custom Python simulations on a high-performance computing cluster. The code is written assuming that the cluster uses the `SLURM` scheduler and that the desired partition name is `general`, so the code may have to be modified to run on your computing resources. It takes no command-line arguments. The code will run 10000 reps of simulations with Q=10, s=0.1 for each of h=0, 0.5 and 1.0 for a population of size N=10000 and a gene conversion rate of 6.55e-05 (10x higher than the base value used for human simulations in the paper; i.e. a 'relative gene conversion rate' of 10).

2. `launchComparisonSlimulations.slurm`: A SLURM script for running the SLiM simulations with the same parameters described immediately above. The output will contain the fraction of sweeps that were pseudo-soft, the mean total frequency of recombinant haplotypes for those pseudo-soft sweeps, and the mean number of gene conversion events per simulation replicate.

3. `compare_to_py.slim`: The SLiM script for running the simulations described immediately above.

4. `combineAllComparisonSimResults.py`: Aggregates all results of the Python simulations into a new directory called `simOutcombined/`. The resulting files contain, among other summaries, the fraction of sweeps that were pseudo-soft, the mean total frequency of recombinant haplotypes for those pseudo-soft sweeps, and the mean number of gene conversion events per simulation replicate.

5. `launchScalingComparisonSlimulations_human.py`: Code for launching human simulations with different values of Q. Uses the `compare_to_py.slim` script to do the actual simulations.

6. `launchScalingComparisonSlimulations_dros.py`: Code for launching Drosophila simulations with different values of Q. Uses the `compare_to_dros.slim` script to do the actual simulations.

7. `compare_to_py_dros.slim`: The SLiM script for running the simulations described immediately above.

8. `commands_for_rescale_comparison.sh`: A bash script that compares the results between the aforementioned SLiM simulations with different values Q.

9. `README.md`: me telling you this.
