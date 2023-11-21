import sys, os

outDir = sys.argv[1]

def getM3FreqFromSamp(haps, positions):
    assert len(haps[0]) < 2 or (len(haps[0]) == 2 and len(positions) == 2)
    if len(haps[0]) == 0:
        return 0
    elif len(haps[0]) == 1:
        if isSweepingAllelePos(positions[0]):
            num = 0 #no m3 found in our sample, only the sweeping allele
        else:
            num = sum([1 for hap in haps if hap[0] == '1']) #sweeping allele is fixed, so count all m3s in sample
        return num
    elif isSweepingAllelePos(positions[1]) and not(isSweepingAllelePos(positions[0])):
        num = 0
        for hap in haps:
            if hap[1] != '1':
                raise ValueError
            if hap[0] == '1':
                num += 1
        return num
    else:
        raise ValueError


def isSweepingAllelePos(pos):
    return pos > 0.4 and pos < 0.6

def getM3PartialFreqFromSamp(haps, positions):
    assert len(positions) == len(haps[0])
    if len(haps[0]) == 2:
        num, denom = 0, 0
        for hap in haps:
            if hap[1] == '1': #beneficial allele
                denom += 1
                if hap[0] == '1':
                    num += 1
    elif len(haps[0]) == 1:
        if not isSweepingAllelePos(positions[0]): # sweeping mutation absent from sample
            num = 0
            denom = 0
        else: # our sweeping mutation present in sample, but m3 was lost
            num = 0
            denom = sum([1 for hap in haps if hap[0] == '1'])
    else:
        num, denom = 0, 0
    return num, denom

def parseSlimResultFile(slimFileName):
    with open(slimFileName, 'rt') as slimFile:
        #print(slimFileName)
        isFixed = 0
        notFixed = 0
        sampMode = 0
        segsites = None
        positions = []
        repNum = int(slimFileName.split("rep_")[1].rstrip(".out"))
        isFixedSoftPop, isFixedSoftSamp, isPartialSoftSamp, restartCount  = 0, 0, 0, 0
        m3freqPop, m3CountSamp, m3PartialCountSamp, m3PartialDenomSamp, n = "NA", "NA", "NA", "NA", "NA"
        sweepStartGen = None
        sojournTime = "NA"
        for line in slimFile:
            if line.startswith("Sweep starts at generation"):
                sweepStartGen = int(line.strip().split()[-1])
            elif line.startswith("NO FIXATION BEFORE SIMULATION END!"):
                notFixed = 1
            elif line.startswith("FIXED at generation"):
                #FIXED at generation 2293
                sweepEndGen = int(line.strip().split()[-1])
                isFixed = 1
            elif line.startswith("LOST - RESTARTING"):
                restartCount += 1
            elif line.startswith("#OUT"):
                out, sampTick, samplCycle, sm, sampledPop, n = line.strip().split() #updated for SLiM 4
                n = int(n)
                sampMode = 1
            elif sampMode == 1:
                sampMode = 2 # ignore "//" line
            elif sampMode == 2:
                segsites = int(line.strip().split()[1])
                if segsites == 0:
                    sampMode = 4
                    haps = []
                else:
                    sampMode = 3
            elif sampMode == 3:
                positions = [float(x) for x in line.strip().split()[1:]]
                sampMode = 4
                haps = []
            elif sampMode == 4:
                haps.append(line.strip())
                if len(haps) == n:
                    sampMode = 5
            elif line.startswith("total m3 frequency:"):
                assert sampMode == 5
                m3freqPop = float(line.strip().split(": ")[-1])
                if m3freqPop > 0:
                    isFixedSoftPop = 1
        #print(isFixed, notFixed)
        if isFixed ^ notFixed:
            goodRep = 1
        else:
            goodRep = 0

        if isFixed:
            assert sweepStartGen != None
            sojournTime = sweepEndGen - sweepStartGen

            m3CountSamp = getM3FreqFromSamp(haps, positions)
            if m3CountSamp > 0:
                isFixedSoftSamp = 1
        elif notFixed:
            m3PartialCountSamp, m3PartialDenomSamp = getM3PartialFreqFromSamp(haps, positions)
            if m3PartialCountSamp > 0:
                isPartialSoftSamp = 1
        if n =="NA":
            goodRep = 0

        return [repNum, goodRep, isFixed, isFixedSoftPop, m3freqPop, isFixedSoftSamp, m3CountSamp, n, isPartialSoftSamp,
                m3PartialCountSamp, m3PartialDenomSamp, sojournTime, restartCount]
    
results = []
goodReps = 0
maxDev = -1
maxDevRep = -1
for fileName in os.listdir(outDir):
    #print(outDir + "/" + fileName)
    result = parseSlimResultFile(outDir + "/" + fileName)
    repNum, goodRep, isFixed, isFixedSoftPop, m3freqPop, isFixedSoftSamp, m3CountSamp, n, isPartialSoftSamp, m3PartialCountSamp, m3PartialDenomSamp, sojournTime, restartCount = result
    goodReps += goodRep
    results.append(result)

    if isFixed:
        sampFreq = m3CountSamp/n
        dev = abs(sampFreq - m3freqPop)
        if dev > maxDev:
            maxDev = dev
            maxDevRep = repNum

for result in results:
    print("\t".join([str(x) for x in result]))

sys.stderr.write(f"{len(results)-goodReps} bad reps out of {len(results)} ({100*(len(results)-goodReps)/len(results)}%)\n")
sys.stderr.write(f"largest deviation between sample and pop m3 frequency: {maxDev}, rep: {maxDevRep}\n")
