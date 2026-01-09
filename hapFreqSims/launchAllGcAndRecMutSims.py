import os
import runCmdAsJob

#modify this function as needed to use your HPC scheduler and resources
def runCmdAsJobWithoutWaitingWithLog(cmd, jobName, launchFile, wallTime, qName, mbMem, logFile):
    """Run the command specified in cmd using the slurm scheduler"""
    with open(launchFile,"w") as f:
        f.write("#!/bin/bash\n")
        f.write("#SBATCH --job-name=%s\n" %(jobName))
        f.write("#SBATCH --time=%s\n" %(wallTime))
        f.write("#SBATCH --partition=%s\n" %(qName))
        f.write("#SBATCH --output=%s\n" %(logFile))
        f.write("#SBATCH --mem=%s\n" %(mbMem))
        f.write("#SBATCH --requeue\n")
        f.write("#SBATCH --export=ALL\n")
        f.write("\n%s\n" %(cmd))
    os.system("sbatch %s" %(launchFile))

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

    jobName="twoLocusSim"
    launchFile="twoLocusSim.slurm"
    wallTime="10-00:00:00"
    qName="general"
    mbMem="8G"

    prefix = f"twoLocusSim_{gcScalar}_{recMutScalar}_{s}_{h}_{Q}"
    if s == 0.001 and h in [0.0, 1.0]:
        repsPerBatch = 1
        numBatches = int(numReps/repsPerBatch)
    elif prefix in intermedPrefixes:
        numBatches = 20
        repsPerBatch = int(numReps/numBatches)
    else:
        numBatches = 1
        repsPerBatch = numReps

    for i in range(numBatches):
        outFile = f"{outDir}/{prefix}_batch_{i}.out"
        cmd = f"python -u simTestRecMut.py {N} {gcScalar} {recMutScalar} {s} {h} {Q} {repsPerBatch} > {outFile}"
        logFile = f"{logDir}/{prefix}_batch_{i}.log"
        runCmdAsJobWithoutWaitingWithLog(cmd, jobName, launchFile, wallTime, qName, mbMem, logFile)
