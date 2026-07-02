#!/bin/bash
#SBATCH --job-name=slimComp
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --output=slim_comparison_results.out
#SBATCH --mem=8G
#SBATCH --requeue
#SBATCH --export=ALL

for d in 0.0 0.5 1.0;
do
    echo $d
    for i in {1..10000} ; do slim -d 'dump_file_name="/tmp/slim_2_locus_test.dump"' -d Q=10 -d popSize=10000 -d unscaledSelCoeff=0.1 -d domCoeff=$d compare_to_py.slim | tail -n 4 ; done > dump.txt ; grep m3 dump.txt | cut -f 4 -d " " | awk '{if ($1 > 0) {freqSum += $1; softies += 1;} tot+=1;} END {print softies/tot, freqSum/softies}'  ; grep recomb dump.txt | cut -f 8 -d " " | aveCols stdin
done
