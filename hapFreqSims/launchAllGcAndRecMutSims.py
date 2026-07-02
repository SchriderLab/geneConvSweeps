import os
import runCmdAsJob

pySimDir = "twoLocusPySims"
mode = "all"

logDir = pySimDir + "/logs"
outDir = pySimDir + "/simOut"

os.system(f"mkdir -p {logDir} {outDir}")

numReps = 100
N=1000000
Q=100
paramCombs = []

for s in [0.001, 0.01, 0.1]:
    for h in [0.0, 0.5, 1.0]:
        for gcScalar in [0, 1, 5, 10, 50]:
            recRateScalar=0
            paramCombs.append((N, gcScalar, recRateScalar, s, h, Q))

        for recRateScalar in [1, 5, 10, 50]:
            gcScalar=0
            paramCombs.append((N, gcScalar, recRateScalar, s, h, Q))


intermedPrefixes = {}
intermedPrefixes['twoLocusSim_0_0_0.01_0.0_100'] = 1
intermedPrefixes['twoLocusSim_1_0_0.01_0.0_100'] = 1
intermedPrefixes['twoLocusSim_5_0_0.01_0.0_100'] = 1
intermedPrefixes['twoLocusSim_10_0_0.01_0.0_100'] = 1
intermedPrefixes['twoLocusSim_0_1_0.01_0.0_100'] = 1
intermedPrefixes['twoLocusSim_0_5_0.01_0.0_100'] = 1
intermedPrefixes['twoLocusSim_0_10_0.01_0.0_100'] = 1
intermedPrefixes['twoLocusSim_0_0_0.01_1.0_100'] = 1


for paramComb in paramCombs:
    N, gcScalar, recMutScalar, s, h, Q = paramComb

    totalGcRate = gcScalar*5e-06
    totalRecMutRate = recMutScalar*1e-08

    jobName="twoLocusSim"
    launchFile="twoLocusSim.slurm"
    wallTime="5-00:00:00"
    qName="general"
    mbMem="8G"

    isInterMed = False
    prefix = f"twoLocusSim_{gcScalar}_{recMutScalar}_{s}_{h}_{Q}"
    if (s == 0.001 and h in [0.0, 1.0]) or (prefix == "twoLocusSim_50_0_0.01_0.0_100"):
        repsPerBatch = 1
        numBatches = int(numReps/repsPerBatch)
    elif prefix in intermedPrefixes:
        numBatches = 20
        repsPerBatch = int(numReps/numBatches)
        isInterMed = True
    else:
        numBatches = 1
        repsPerBatch = numReps

    for i in range(numBatches):
        outFile = f"{outDir}/{prefix}_batch_{i}.out"
        cmd = f"python -u simTestRecMut.py {N} {totalGcRate} {totalRecMutRate} {s} {h} {Q} {repsPerBatch} > {outFile}"
        logFile = f"{logDir}/{prefix}_batch_{i}.log"
        if mode == "all" or (mode == "interMedOnly" and isInterMed):
            runCmdAsJob.runCmdAsJobWithoutWaitingWithLog(cmd, jobName, launchFile, wallTime, qName, mbMem, logFile)
