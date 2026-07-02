echo "comparison of fraction of sweeps that are pseudo-soft -- human"
for s in 0.01 0.1 ; do echo "s=$s" ; for Q in 1 10 ; do echo "Q=$Q" ; for h in 0.0 0.5 1.0 ; do echo "h=$h" ; grep "total m3" slimSimsForRescsalingTest/PopSize_10000/simOut/constPop_10000_10_${s}_${h}_${Q}_batch_* | cut -f 4 -d " " | awk '{if ($1 > 0) softies += 1; tot += 1} END {print softies/tot}' ; done ; echo "" ; done ; printf "\n" ; done

echo "comparison of recombinant haplotype frequency -- human"
for s in 0.01 0.1 ; do echo "s=$s" ; for Q in 1 10 ; do echo "Q=$Q" ; for h in 0.0 0.5 1.0 ; do echo "h=$h" ; grep "total m3" slimSimsForRescsalingTest/PopSize_10000/simOut/constPop_10000_10_${s}_${h}_${Q}_batch_* | cut -f 4 -d " " | awk '{if ($1 > 0) {softies += 1; freq += $1}} END {print freq/softies}' ; done ; echo "" ; done ; printf "\n" ; done

echo "comparison of fraction of sweeps that are pseudo-soft -- dmel"
for s in 0.01 0.1 ; do echo "s=$s" ; for Q in 10 100 ; do echo "Q=$Q" ; for h in 0.0 0.5 1.0 ; do echo "h=$h" ; grep "total m3" slimSimsForRescsalingTest_dros/PopSize_1720600/simOut/constPop_1720600_1_${s}_${h}_${Q}_batch_* | cut -f 4 -d " " | awk '{if ($1 > 0) softies += 1; tot += 1} END {print softies/tot}' ; done ; echo "" ; done ; printf "\n" ; done

echo "comparison of recombinant haplotype frequency -- dmel"
for s in 0.01 0.1 ; do echo "s=$s" ; for Q in 10 100 ; do echo "Q=$Q" ; for h in 0.0 0.5 1.0 ; do echo "h=$h" ; grep "total m3" slimSimsForRescsalingTest_dros/PopSize_1720600/simOut/constPop_1720600_1_${s}_${h}_${Q}_batch_* | cut -f 4 -d " " | awk '{if ($1 > 0) {softies += 1; freq += $1}} END {print freq/softies}' ; done ; echo "" ; done ; printf "\n" ; done
