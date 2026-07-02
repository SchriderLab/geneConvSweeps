import os
sys.path.insert(0, "../hapFreqSims")
import runCmdAsJob

pySimDir = "twoLocusPySimsForComparison"

logDir = pySimDir + "/logs"
outDir = pySimDir + "/simOut"

os.system(f"mkdir -p {logDir} {outDir}")

numReps = 10000
N=10000
Q=10
paramCombs = []

for s in [0.1]:
    for h in [0.0, 0.5, 1.0]:
        for gcScalar in [10]:
            recRateScalar=0
            paramCombs.append((N, gcScalar, recRateScalar, s, h, Q))

gcRate=6.55e-06
recMutRate=1e-8
for paramComb in paramCombs:
    N, gcScalar, recMutScalar, s, h, Q = paramComb

    jobName="2LocCompSim"
    launchFile="2LocCompSim.slurm"
    wallTime="5-00:00:00"
    qName="general"
    mbMem="8G"

    prefix = f"twoLocusSim_{gcScalar}_{recMutScalar}_{s}_{h}_{Q}"
    repsPerBatch = 100
    numBatches = int(numReps/repsPerBatch)

    for i in range(numBatches):
        outFile = f"{outDir}/{prefix}_batch_{i}.out"
        cmd = f"python -u ../hapFreqSims/simTestRecMut.py {N} {gcScalar*gcRate} {recMutScalar*recMutRate} {s} {h} {Q} {repsPerBatch} > {outFile}"
        logFile = f"{logDir}/{prefix}_batch_{i}.log"
        runCmdAsJob.runCmdAsJobWithoutWaitingWithLog(cmd, jobName, launchFile, wallTime, qName, mbMem, logFile)
