import os
sys.path.insert(0, "../hapFreqSims")
import runCmdAsJob


N=1720600
simDir = f"slimSimsForRescsalingTest_dros/PopSize_{N}"

logDir = simDir + "/logs"
outDir = simDir + "/simOut"
dumpDir = simDir + "/dump"

os.system(f"mkdir -p {logDir} {outDir} {dumpDir}")

numReps = 1000

paramCombs = []
for Q in [10, 100]:
    for s in [0.1, 0.01]:
        for h in [0.0, 0.5, 1.0]:
            for gcScalar in [1]:  # this is hard-coded in the SLiM script anyway
                assert gcScalar == 1
                recRateScalar=0
                paramCombs.append((N, gcScalar, s, h, Q))



for paramComb in paramCombs:
    N, gcScalar, s, h, Q = paramComb

    jobName="SlimRescale_dros"
    launchFile="SlimRescale_dros.slurm"
    wallTime="10-00:00:00"
    qName="general"
    mbMem="8G"

    prefix = f"constPop_{N}_{gcScalar}_{s}_{h}_{Q}"
    repsPerBatch = 100
    numBatches = int(numReps/repsPerBatch)

    for i in range(numBatches):
        cmd = "eval \"$(conda shell.bash hook)\"\nconda activate stdpopsim\n"
        cmd += f"for j in {{0..{repsPerBatch-1}}} ; do slim -d \"dump_file_name='{dumpDir}/{prefix}_batch_{i}_rep_$j.dump'\" -d Q={Q} -d popSize={N} -d unscaledSelCoeff={s} -d domCoeff={h} compare_to_py_dros.slim > {outDir}/{prefix}_batch_{i}_rep_$j.out ; done"
        logFile = f"{logDir}/{prefix}_batch_{i}.log"
        runCmdAsJob.runCmdAsJobWithoutWaitingWithLog(cmd, jobName, launchFile, wallTime, qName, mbMem, logFile)
