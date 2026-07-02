import sys, os

def runCmdAsJobWithoutWaitingWithLog(cmd, jobName, launchFile, wallTime, qName, mbMem, logFile):
    with open(launchFile,"w") as f:
        f.write("#!/bin/bash\n")
        f.write("#SBATCH --job-name=%s\n" %(jobName))
        f.write("#SBATCH --time=%s\n" %(wallTime))
        f.write("#SBATCH --partition=%s\n" %(qName))
        f.write("#SBATCH --output=%s\n" %(logFile))
        f.write("#SBATCH --mem=%s\n" %(mbMem))
        f.write("#SBATCH --requeue\n")
        f.write("#SBATCH --export=ALL\n")
        if qName == "volta-gpu":
            f.write("#SBATCH --gres=gpu:1\n#SBATCH --qos=gpu_access\n\n")
            f.write("unset OMP_NUM_THREADS\nSIMG_PATH=/nas/longleaf/apps/tensorflow_py3/1.12.0/simg\n")
            f.write("SIMG_NAME=tensorflow1.12.0-py3-cuda9.0-ubuntu16.04.simg\n")
            f.write("DATA_PATH={}\n".format(os.getcwd()))
            f.write("\nsingularity exec --nv -B /pine -B /proj $SIMG_PATH/$SIMG_NAME bash -c \"%s\"\n" %(cmd))
        else:
            f.write("\n%s\n" %(cmd))
    os.system("sbatch %s" %(launchFile))

def main():
    cmd, jobName, launchFile, wallTime, qName, mbMem, logFile = sys.argv[1:]
    runCmdAsJobWithoutWaitingWithLog(cmd, jobName, launchFile, wallTime, qName, mbMem, logFile)

if __name__ == "__main__":
    main()
