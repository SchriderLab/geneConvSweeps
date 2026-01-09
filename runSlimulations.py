import sys, os

species = sys.argv[1]

scriptDir = "slimScriptsInjected"
logDir = "slimLogs"
outDir = "slimOutput"
partitionName = "general" #replace if you want to run the slurm job on a different partition


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


os.system(f"mkdir -p {outDir} {logDir}")
repRange = "0-9" #we are going to launch 10 jobs

for scriptFileName in [x for x in os.listdir(scriptDir) if species in x]:
    scriptFileNoExt = scriptFileName.rstrip(".slim")
    outPath = outDir + "/" + scriptFileNoExt
    logPath = logDir + "/" + scriptFileNoExt
    os.system(f"mkdir -p {outPath} {logPath}")
    logPath += "/repBatch_%a.log"
    jobName = "gcSweep"
    launchFileName = "gcSweep.slurm"

    #each job runs 100 replicates
    cmd = f"""
repStart=$(( SLURM_ARRAY_TASK_ID*100 ))
repEnd=$(( repStart + 99 ))
for (( i=$repStart; i<=$repEnd; i++ ))
do
    echo "starting rep $i"
    slim -d rep_number=$i {scriptDir}/{scriptFileName} > {outPath}/rep_$i.out
done"""

    if species == "DroMel":
        timeLimit="5-00:00:00"
    else:
        timeLimit="12:00:00"
    runCmdAsJobArrayWithoutWaitingWithLog(cmd, jobName, launchFileName, timeLimit, partitionName, "8G", logPath, repRange)
