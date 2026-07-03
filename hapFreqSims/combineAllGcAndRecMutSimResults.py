import os
import sys
import math
import ast
import numpy as np

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



def parseHapCounterGetHapFreqs(line):
    line = line.strip()
    line = line.lstrip("Counter(")
    line = line.rstrip(")")
    c = ast.literal_eval(line)
    vals = list(c.values())
    denom = sum(vals)
    freqs = [x/denom for x in vals]
    return freqs



def getHapStats(hapFreqs):
    hapFreqs.sort(reverse=True)
    maf = 1-(hapFreqs[0])
    nHaps = len(hapFreqs)

    if nHaps > 1:
        squareFreqs = [x**2 for x in hapFreqs]
        h1 = sum(squareFreqs)
        h2 = sum(squareFreqs[1:])
        h2h1 = h2/h1
        #print("hapFreqs, squareFreqs, h1, h2, h2h1, maf, nHaps")
        #print(hapFreqs, squareFreqs, h1, h2, h2h1, maf, nHaps)
    else:
        h1 = 1.0
        h2 = 0.0
        h2h1 = 0.0

    return maf, h2h1, nHaps



def getAvgHapStats(allHapFreqs):
    mafs, h2h1s, hapcounts = [], [], []
    for hapFreqs in allHapFreqs:
        maf, h2h1, hapcount = getHapStats(hapFreqs)
        if hapcount > 1:
            mafs.append(maf)
            h2h1s.append(h2h1)
        hapcounts.append(hapcount)

    #print(h2h1s)
    return np.mean(mafs), np.mean(h2h1s), np.mean(hapcounts)



def readResFile(filePath):
    nReps = 0
    sampMode = 0
    batchHapFreqsPop, batchHapFreqsSamp = [], []
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
            elif line.startswith("avg number of gene conversion events per replicate"):
                geneConvEvents = float(line.strip().split()[-1])
            elif line.startswith("avg pop marker mut freq in soft sweeps"):
                softMarkFreqPop, softMarkFreqSamp = parseResSummaryLine(line)
            elif line.startswith("Counter("):
                if sampMode == 0:
                    hapFreqsPop = parseHapCounterGetHapFreqs(line)
                    batchHapFreqsPop.append(hapFreqsPop)
                elif sampMode == 1:
                    hapFreqsSamp = parseHapCounterGetHapFreqs(line)
                    batchHapFreqsSamp.append(hapFreqsSamp)
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
        softKSamp,
        softMarkFreqPop,
        softMarkFreqSamp,
        geneConvEvents,
        batchHapFreqsPop,
        batchHapFreqsSamp
    )
    return retVals



def combineResultsFilesWithPrefix(outDir, prefix):
    fnames = getFileNamesWithPrefix(outDir, prefix)
    resTots = [0]*11
    denoms = [0]*11
    softOnlyIndices = [3, 4, 5, 6, 9, 10]

    allHapFreqsPop, allHapFreqsSamp = [], []
    for fname in fnames:
        #print(fname)
        skip = False
        try:
            res = readResFile(outDir + "/" + fname)
            allHapFreqsPop.extend(res[-2])
            allHapFreqsSamp.extend(res[-1])
        except Exception as e:
            skip = True
            sys.stderr.write(f"skipping {fname} which is incomplete\n")

        if not skip:
            nReps = res[0]
            for i in range(1, len(res)-2):
                if i in softOnlyIndices:
                    if not math.isnan(res[i]):
                        softFrac = res[2 - (i % 2)]
                        nSoft = softFrac*nReps
                        denoms[i-1] += nSoft
                        resTots[i-1] += res[i]*nSoft
                else:
                    resTots[i-1] += res[i]*nReps
                    denoms[i-1] += nReps

    for i in range(len(resTots)):
        if denoms[i] == 0:
            resTots[i] = 'NA'
        else:
            resTots[i] /= denoms[i]

    avgMafPop, avgH2H1Pop, avgNHapsPop = getAvgHapStats(allHapFreqsPop)
    avgMafSamp, avgH2H1Samp, avgNHapsSamp = getAvgHapStats(allHapFreqsSamp)

    outLines=[]
    outLines.append(f"fraction of soft sweeps in pop: {resTots[0]}; fraction of soft sweeps in sample: {resTots[1]}")

    # found a bug for the averaging of these things within the simulation batches, but found it is easier to just re-calculate them here
    #outLines.append(f"avg pop h2/h1 of soft sweeps: {resTots[2]}; avg sample h2/h1 of soft sweeps: {resTots[3]}")
    #outLines.append(f"avg pop minor hap freq of soft sweeps: {resTots[4]}; avg sample minor hap freq of soft sweeps: {resTots[5]}")
    #outLines.append(f"avg number of participating haps in pop: {resTots[6]}; avg number of participating haps in samp: {resTots[7]}")

    outLines.append(f"avg pop h2/h1 of soft sweeps: {avgH2H1Pop}; avg sample h2/h1 of soft sweeps: {avgH2H1Samp}")
    outLines.append(f"avg pop minor hap freq of soft sweeps: {avgMafPop}; avg sample minor hap freq of soft sweeps: {avgMafSamp}")
    outLines.append(f"avg number of participating haps in pop: {avgNHapsPop}; avg number of participating haps in samp: {avgNHapsSamp}")
    outLines.append(f"avg pop marker mut freq in soft sweeps: {resTots[8]}; avg samp marker mut freq in soft sweeps: {resTots[9]}")
    outLines.append(f"avg number of gene conversion events per replicate: {resTots[10]}")

    return outLines




for paramComb in paramCombs:
    N, gcScalar, recMutScalar, s, h, Q = paramComb

    prefix = f"twoLocusSim_{gcScalar}_{recMutScalar}_{s}_{h}_{Q}"
    #if prefix != "twoLocusSim_0_1_0.01_0.0_100":
    #    continue

    outFilePath = combinedDir + "/" + prefix + ".out"

    outLines = combineResultsFilesWithPrefix(outDir, prefix)

    sys.stderr.write(f"combining results for {outFilePath}\n")
    with open(outFilePath, "w") as f:
        for l in outLines:
            f.write(l + "\n")
