import os
import sys
import math
import ast

pySimDir = "twoLocusPySims"

outDir = pySimDir + "/simOut"
combinedDir = pySimDir + "/simOutCombined"

os.system(f"mkdir -p {combinedDir}")

numReps = 100
N=1000000
Q=100
paramCombs = []

for s in [0.001, 0.01, 0.1]:
    for h in [0.0, 0.5, 1.0]:
        for gcScalar in [0, 1, 5, 10, 50]:
            recRateScalar=0
            paramCombs.append((N, gcScalar, recRateScalar, s, h, Q))

        #for recRateScalar in [1, 10, 100, 1000]:
        for recRateScalar in [1, 5, 10, 50]:
            gcScalar=0
            paramCombs.append((N, gcScalar, recRateScalar, s, h, Q))

def getFileNamesWithPrefix(outDir, prefix):
    fnames = []
    for fname in os.listdir(outDir):
        if fname.split("_batch")[0] == prefix:
            fnames.append(fname)
    return fnames

def parseResSummaryLine(resSummaryLine):
    resSummaryLine = resSummaryLine.rstrip()
    n1 = float(resSummaryLine.split(":")[1].split(";")[0])
    n2 = float(resSummaryLine.split(" ")[-1])
    return n1, n2

def parseHapCounterGetMajHapFreq(line):
    line = line.strip()
    line = line.lstrip("Counter(")
    line = line.rstrip(")")
    c = ast.literal_eval(line)
    vals = list(sorted(c.values(), reverse=True))
    denom = sum(vals)
    majorHapFreq = max(vals)/denom
    return majorHapFreq

def readResFile(filePath):
    nReps = 0
    sampMode = 0
    with open(filePath) as f:
        for line in f:
            if line.startswith("fraction of soft sweeps in pop"):
                softFracPop, softFracSamp = parseResSummaryLine(line)
            elif line.startswith("avg pop h2/h1 of soft sweeps"):
                softH21Pop, softH21Samp = parseResSummaryLine(line)
            elif line.startswith("avg minor hap pop freq of soft sweeps"):
                softFreqPop, softFreqSamp = parseResSummaryLine(line)
            elif line.startswith("avg number of participating haps in pop"):
                softKPop, softKSamp = parseResSummaryLine(line)
            elif line.startswith("Counter("):
                if sampMode == 0:
                    softFreqPop = 1 - parseHapCounterGetMajHapFreq(line)
                elif sampMode == 1:
                    softFreqSamp = 1 - parseHapCounterGetMajHapFreq(line)
                else:
                    raise Exception(f"more than two counter lines encountered for rep {nReps} in {filePath}")
                sampMode += 1
            elif line.startswith("done with rep"):
                nReps += 1
                sampMode = 0
    retVals = (
        nReps,
        softFracPop,
        softFracSamp,
        softH21Pop,
        softH21Samp,
        softFreqPop,
        softFreqSamp,
        softKPop,
        softKSamp
    )
    return retVals

def combineResultsFilesWithPrefix(outDir, prefix):
    fnames = getFileNamesWithPrefix(outDir, prefix)
    totReps = 0
    resTots = [0]*8
    for fname in fnames:
        res = readResFile(outDir + "/" + fname)
        nReps = res[0]
        for i in range(1, len(res)):
            if i in [3, 4, 5, 6]:
                if not math.isnan(res[i]):
                    softFrac = res[2 - (i % 2)]
                    resTots[i-1] += res[i]*softFrac*nReps
            else:
                resTots[i-1] += res[i]*nReps
        totReps += nReps
    for i in range(len(resTots)):
        resTots[i] /= totReps
    outLines=[]
    outLines.append(f"fraction of soft sweeps in pop: {resTots[0]}; fraction of soft sweeps in sample: {resTots[1]}")
    outLines.append(f"avg pop h2/h1 of soft sweeps: {resTots[2]}; avg sample h2/h1 of soft sweeps: {resTots[3]}")
    outLines.append(f"avg pop minor hap freq of soft sweeps: {resTots[4]}; avg sample minor hap freq of soft sweeps: {resTots[5]}")
    outLines.append(f"avg number of participating haps in pop: {resTots[6]}; avg number of participating haps in samp: {resTots[7]}")

    return outLines

for paramComb in paramCombs:
    N, gcScalar, recMutScalar, s, h, Q = paramComb

    prefix = f"twoLocusSim_{gcScalar}_{recMutScalar}_{s}_{h}_{Q}"
    outFilePath = combinedDir + "/" + prefix + ".out"

    outLines = combineResultsFilesWithPrefix(outDir, prefix)

    sys.stderr.write(f"combining results for {outFilePath}\n")
    with open(outFilePath, "w") as f:
        for l in outLines:
            f.write(l + "\n")
