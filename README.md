# geneConvSweeps
This repository contains code for simulation study exaiming the impact of gene conversion on the "softness" of selective sweeps, used for the analyses presented in a forthcoming preprint.

# Simulation pipeline contents
This repository contains four python scripts that generate, modify, and run SLiM scripts (for use with SLiM version 4.0.1) simulating a two-locus hitchhiking model with gene conversion events allowed at the selected locus but not the linked locus. This pipeline requires that `stdpopsim` version 0.2.0 and all of its dependencies be installed. (This is pretty easy using `conda` following the instructions on https://popsim-consortium.github.io/stdpopsim-docs/stable/installation.html). These scripts write all of their output to directories that they will create within the current working directory, so if you want to write them somewhere else you will have to modify them slightly.

1. `buildAllSlimScripts.py`: This script uses `stdpopsim` to generate all of the SLiM scripts which we will then modify to run our sweep simulations under the desired demographic models (using some *Arabidopsis*, *Drosophila*, and human models from the `stdpopsim` catalog).

2. `injectAllSlimScripts.py`: This code modifies the SliM scripts generated by `stdpopsim` via the above script. These SLiM scripts are modified to contain selective sweeps with gene conversion occurring only at the selected locus. Note that this step could be simplified substantially by taking advantage of `stdpopsim`'s ability to condition on sweeps occuring at a specified time (by adding code for this to `buildAllSlimScripts.py`), but the development of this pipeline began prior to the incorporation of that functionality into `stdpopsim`.

3. `runSlimulations.py`: This code actually launches the simulation jobs to a high-performance computing cluster. The code is written assuming that the cluster uses the `SLURM` scheduler and that the desired partition name is `general`, so the code may have to be modified to run on your computing resources. It can be run for each species as follows:
   `python runSlimulations.py HomSap`
   `python runSlimulations.py AraTha`
   `python runSlimulations.py DroMel`

4. `parseOutputForSpecies.py`: This code parses the output from each SLiM simulation and writes information about the simulation's outcomes into a tabular format that can be read by the analysis notebook described below. This can be run for each species as follows:
   `python parseOutputForSpecies.py HomSap`
   `python parseOutputForSpecies.py DroMel`
   `python parseOutputForSpecies.py AraTha`

# Notebook with downstream analysis
Need to add this next and add a brief description of it.
