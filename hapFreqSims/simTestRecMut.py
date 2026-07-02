import sys
import random, copy, collections
import numpy as np

class Population:

    def selectParent(self):
        return random.choices(self.genomes, weights=self.fitnesses, k=1)[0]
    
    def selectRandomGenomeIndex(self):
        return random.randrange(len(self.genomes))
        
    def __init__(self, genomes=[], N=0, s=0, h=0.5):
        if (genomes and N) or not (genomes or N):
            raise ValueError('Must define genomes list or N, not both')
        if N > 0:
            genomes = []
            for i in range(N):
                genomes.append(Genome(Chromosome(self), Chromosome(self)))

        self.genomes = genomes
        self.fitnesses = [1]*len(genomes)
        self.s = s
        self.h = h
        self.currSelMutCount = 0
           
    def updateFitnesses(self):
        for i in range(len(self.genomes)):
            self.fitnesses[i] = self.genomes[i].getFitness(self)
            
    def updateFitnessOfIndexedGenome(self, i):
        self.fitnesses[i] = self.genomes[i].getFitness(self)
        
    def getSelMutFreq(self):
        num, denom = 0, 0
        for genome in self.genomes:
            num += genome.chromosome1.hasSelMut() + genome.chromosome2.hasSelMut()
            denom += 2
        return num/denom
    
    def getMarkerMutFreq(self, n=None, maf=True):
        num, denom = 0, 0
        if n:
            genomes = random.sample(self.genomes, n)
            print(f"sample genotypes for {n} individuals: {n}")
            genotypes = []
        else:
            genomes = self.genomes
        for genome in genomes:
            num += genome.chromosome1.hasMarkerMut() + genome.chromosome2.hasMarkerMut()
            denom += 2
            if n:
                genotypes.append((genome.chromosome1.markerMut, genome.chromosome1.selMut))
                genotypes.append((genome.chromosome2.markerMut, genome.chromosome2.selMut))
        if n:
            c = collections.Counter(genotypes)
            print("marker geno counter:", c)
            print(len(c))
        freq = num/denom
        if maf and freq > 0.5:
            freq = 1-freq

        return freq
    
    def countParticipatingHaps(self, n=None):
        if n:
            genomes = random.sample(self.genomes, n)
        else:
            genomes = self.genomes
        selAlleles = []
        for genome in genomes:
            selAlleles.extend([genome.chromosome1.selMut, genome.chromosome2.selMut])
        c = collections.Counter(selAlleles)
        print(c)
        print(len(c))
        return len(c)

    def getH2OverH1(self, n=None):
        if n:
            genomes = random.sample(self.genomes, n)
        else:
            genomes = self.genomes
        selAlleles = []
        for genome in genomes:
            selAlleles.extend([genome.chromosome1.selMut, genome.chromosome2.selMut])
        c = collections.Counter(selAlleles)

        hapFreqs = []
        numChroms=2*len(genomes)
        for hapID in c:
            hapCount = c[hapID]
            hapFreq = hapCount/numChroms
            hapFreqs.append(hapFreq)
        
        if len(hapFreqs) > 1:
            hapFreqs.sort(reverse=True)
            #print(hapFreqs, sum(hapFreqs))
            h1 = sum([x**2 for x in hapFreqs])
            h12 = (hapFreqs[0]+hapFreqs[1])**2 + sum([x**2 for x in hapFreqs[2:]])
            h2 = sum([x**2 for x in hapFreqs[1:]])
            #print(h1, h2)
            mhf = 1-hapFreqs[0]
        else:
            h1, h2 = 1.0, 0.0
            mhf = 0.0
        return h2/h1, mhf



class Chromosome:
    
    def addMutation(self):
        if not self.hasSelMut(): #not allowing chroms with B allele to mutate to B allele
            self.population.currSelMutCount += 1
            self.selMut = self.population.currSelMutCount
    
    def removeMarkerMut(self):
        self.markerMut = 0
            
    def __init__(self, pop, selMut=0, markerMut=1):
        self.selMut = selMut
        self.markerMut = markerMut
        self.population = pop
        
    def hasSelMut(self):
        if self.selMut == 0:
            return 0
        else:
            return 1
    
    def hasMarkerMut(self):
        return self.markerMut
        


class Genome:
    def getFitness(self, pop):
        mutCopies = self.chromosome1.hasSelMut() + self.chromosome2.hasSelMut()
        if mutCopies == 2:
            fitness = 1 + pop.s
        elif mutCopies == 1:
            fitness = 1 + pop.h*pop.s
        else:
            fitness = 1
        return fitness
    
    def addSelMutAndRemoveMarker(self):
        chrom = random.choice([self.chromosome1, self.chromosome2])
        chrom.addMutation()
        chrom.removeMarkerMut()
        
    def __init__(self, chromosome1, chromosome2):
        self.chromosome1 = chromosome1
        self.chromosome2 = chromosome2
        
    def clone(self):
        newChrom1 = Chromosome(self.chromosome1.population, selMut=self.chromosome1.selMut, markerMut=self.chromosome1.markerMut)
        newChrom2 = Chromosome(self.chromosome2.population, selMut=self.chromosome2.selMut, markerMut=self.chromosome2.markerMut)
        return Genome(newChrom1, newChrom2)



def fmtMarker(m):
    assert m in [0, 1]
    if m == 1:
        return "a"
    # initially the selected mutation occurs on a haplotype without the maker mutation in this sim
    # we call the marker-less allele the A allele in the paper
    elif m == 0:
        return "A"



def fmtSelMut(m):
    assert m in [0, 1]
    if m == 1:
        return "B"  # our selected allele
    elif m == 0:
        return "b"  # ancestral allele at the selected site



def doConversionAtSelSite(donorChrom, recipChrom):
    conversionType = f"{fmtMarker(donorChrom.hasMarkerMut())}{fmtSelMut(donorChrom.hasSelMut())}->{fmtMarker(recipChrom.hasMarkerMut())}{fmtSelMut(recipChrom.hasSelMut())}"
    if donorChrom.selMut > 0 and recipChrom.selMut == 0:
        recipChrom.population.currSelMutCount += 1
        recipChrom.selMut = recipChrom.population.currSelMutCount
    elif donorChrom.selMut == 0:
        recipChrom.selMut = 0
    return conversionType



def makeRecombinantGametesForParent(parent, totalGcRate):
    newParent = parent.clone()
    chrom1Converted = None
    chrom2Converted = None
    gcEvent = 0
    # doubling totalGcRate here because our rate should signify the probability that the *chosen gamete* is converted
    # a more efficient way would be to just generate a single (possibly recombinant) gamete rather than generate both and
    # randomly pick one later in makeRecombinantOffspringChroms, but this is what i did...
    if random.random() < 2*totalGcRate:
        if random.randint(0, 1):
            conversionType = doConversionAtSelSite(newParent.chromosome2, newParent.chromosome1)
            chrom1Converted = conversionType
        else:
            conversionType = doConversionAtSelSite(newParent.chromosome1, newParent.chromosome2)
            chrom2Converted = conversionType
    
    return (newParent.chromosome1, chrom1Converted), (newParent.chromosome2, chrom2Converted)



def recordConversionTypes(conv1, conv2):
    convs = {}
    for c in conv1, conv2:
        if c != None:
            if not c in convs:
                convs[c] = 0
            convs[c] += 1
    return convs



def addGcEvents(totalGcEvents, newGcEvents):
    for key in newGcEvents:
        if key in totalGcEvents:
            totalGcEvents[key] += newGcEvents[key]
        else:
            totalGcEvents[key] = newGcEvents[key]



def countTotalGcEvents(gcEvents):
    tot = 0
    convToNew, convToOrig, miscConv = [], [], []
    for k in gcEvents:
        tot += gcEvents[k]

    #    # for checking rate of origination of old and new haps via gene conv:
    #    "aB->Ab"
    #    if k[1] == "B" and k[4:] == "ab":
    #        convToNew.append((gcEvents[k], k))
    #    elif k[1] == "B" and k[4:] == "Ab":
    #        convToOrig.append((gcEvents[k], k))
    #    else:
    #        miscConv.append((gcEvents[k], k))

    #print("gene conversion events that yield a new hap")
    #for v, k in reversed(sorted(convToNew)):
    #    print(k, v)
    #print("gene conversion events that yield the orig hap")
    #for v, k in reversed(sorted(convToOrig)):
    #    print(k, v)
    #print("misc other gene conversion events")
    #for v, k in reversed(sorted(miscConv)):
    #    print(k, v)
        
    return tot



def makeRecombinantOffspringChroms(parent1, parent2, totalGcRate):
    parent1Gametes = makeRecombinantGametesForParent(parent1, totalGcRate)
    parent1Gamete, isConverted1 = random.choice(parent1Gametes)
    
    parent2Gametes = makeRecombinantGametesForParent(parent2, totalGcRate)
    parent2Gamete, isConverted2 = random.choice(parent2Gametes)
    
    return parent1Gamete, parent2Gamete, recordConversionTypes(isConverted1, isConverted2)
        


def makeOffspringGenome(parent1, parent2, totalGcRate, recMutRate):
    offspringChrom1, offspringChrom2, gcEvents = makeRecombinantOffspringChroms(parent1, parent2, totalGcRate)
    
    # possibility of recurrent mutation to the adaptive allele
    if random.random() < recMutRate:
        offspringChrom1.addMutation()
    if random.random() < recMutRate:
        offspringChrom2.addMutation()
        
    offspringGenome = Genome(offspringChrom1, offspringChrom2)
    return offspringGenome, gcEvents



def runSimUntilFixation(Ngen, Nt, pop, totalGcRate, recMutRate):
    done = 0
    popH21, sampleH21 = None, None
    established, nearFixation = False, False
    totalGcEvents = {}  # counts total GC events during the simulation, including during failed sweeps that were restarted
    while not done:
        for gen in range(Ngen):
            currN = Nt[gen]
            establishmentThreshold = (1/pop.s)/currN
            #establishmentThreshold = 0.05
            nearFixationThreshold = 1-establishmentThreshold
            genomes = []

            for offspring in range(currN):
                parent1 = pop.selectParent()
                parent2 = pop.selectParent()
                offspringGenome, gcEvents = makeOffspringGenome(parent1, parent2, totalGcRate, recMutRate)
                genomes.append(offspringGenome)
                addGcEvents(totalGcEvents, gcEvents)

            pop = Population(genomes=genomes, s=pop.s, h=pop.h)
            pop.updateFitnesses()
            selMutFreq = pop.getSelMutFreq()
            if selMutFreq == 0:
                #print(f"selected mutation lost at generation: {gen}. restarting")
                #init again:
                established = False
                nearFixation = False
                pop = Population(N=Nanc, s=pop.s, h=pop.h)
                targetGenomeIndex = pop.selectRandomGenomeIndex()
                pop.genomes[targetGenomeIndex].addSelMutAndRemoveMarker()
                pop.updateFitnessOfIndexedGenome(targetGenomeIndex)
                break
            elif (not established) and selMutFreq > establishmentThreshold:
                if verbose:
                    print(f"selected mutation established (freq={selMutFreq}) at generation: {gen}")
                established=True
            elif (not nearFixation) and selMutFreq > nearFixationThreshold:
                if verbose:
                    print(f"selected mutation near fixation (freq={selMutFreq}) at generation: {gen}")
                nearFixation=True
            if selMutFreq == 1:
                if verbose:
                    print(f"selected mutation fixed at generation: {gen}")
                popMarkerFreq, sampleMarkerFreq = pop.getMarkerMutFreq(maf=False), pop.getMarkerMutFreq(n=100, maf=False)
                popH21, popMhf = pop.getH2OverH1()
                sampleH21, sampleMhf = pop.getH2OverH1(n=100)
                popNumParticipatingHaps, sampleNumParticipatingHaps = pop.countParticipatingHaps(), pop.countParticipatingHaps(n=100)
                if (sampleNumParticipatingHaps) > 1:
                    print("soft in samp")
                if (popNumParticipatingHaps) > 1:
                    print("soft in pop")
                done = 1
                break
    assert popH21 != None and sampleH21 != None
    return gen, popMarkerFreq, sampleMarkerFreq, popH21, popMhf, sampleH21, sampleMhf, popNumParticipatingHaps, sampleNumParticipatingHaps, totalGcEvents

def runReplicate(Nanc, Ngen, totalGcRate, recMutRate, s, h):

    #initialize the population here
    #only constant-size pops currently supported
    Nt = [Nanc]*Ngen
    pop = Population(N=Nanc, s=s, h=h)

    #then randomly pick a genome
    targetGenomeIndex = pop.selectRandomGenomeIndex()

    #then add the mutation there, give it the marker mutation too
    pop.genomes[targetGenomeIndex].addSelMutAndRemoveMarker()

    #then update its fitness
    pop.updateFitnessOfIndexedGenome(targetGenomeIndex)

    return runSimUntilFixation(Ngen, Nt, pop, totalGcRate, recMutRate)

baseNanc, baseGcRate, baseRecMutRate, s, h, Q, numReps = sys.argv[1:]
verbose=True

baseGcRate = float(baseGcRate)
baseRecMutRate = float(baseRecMutRate)
s=float(s)
h=float(h)
Q=float(Q)
numReps = int(numReps)

baseNanc=float(baseNanc)
Nanc = int(baseNanc/Q)
Ngen = Nanc*10
totalGcRate = baseGcRate*Q
recMutRate = baseRecMutRate*Q
sys.stderr.write(f"running {numReps} reps with N={baseNanc}, s={s}, h={h}, 4Ng={4*Nanc*totalGcRate} and 4Nu={4*Nanc*recMutRate}\n")
s = s*Q


gens, popMarkerFreqs, sampMarkerFreqs, popH21s, popMhfs, sampleH21s, sampleMhfs, popHaps, sampleHaps, gcEventCounts = [], [], [], [], [], [], [], [], [], []

for rep in range(numReps):
    print(f"starting rep {rep}")
    gen, popMarkerFreq, sampMarkerFreq, popH21, popMhf, sampleH21, sampleMhf, popNumParticipatingHaps, sampleNumParticipatingHaps, totalGcEvents = runReplicate(
        Nanc, Ngen, totalGcRate, recMutRate, s, h)
    print(f"done with rep {rep}")

    gens.append(gen)
    if popH21 != 0:
        popH21s.append(popH21)
        popMhfs.append(popMhf)
        popMarkerFreqs.append(popMarkerFreq)
    if sampleH21 != 0:
        sampleH21s.append(sampleH21)
        sampleMhfs.append(sampleMhf)
        sampMarkerFreqs.append(sampMarkerFreq)
    popHaps.append(popNumParticipatingHaps)
    sampleHaps.append(sampleNumParticipatingHaps)
    gcEventCounts.append(countTotalGcEvents(totalGcEvents))

    print("")

popSofts = [x for x in popHaps if x > 1]
sampleSofts = [x for x in sampleHaps if x > 1]

print(f"fraction of soft sweeps in pop: {len(popSofts)/len(popHaps)}; fraction of soft sweeps in sample: {len(sampleSofts)/len(sampleHaps)}")
print(f"avg pop h2/h1 of soft sweeps: {np.mean(popH21s)}; avg sample h2/h1 of soft sweeps: {np.mean(sampleH21s)}")
print(f"avg minor hap pop freq of soft sweeps: {np.mean(popMhfs)}; avg minor hap sample freq of soft sweeps: {np.mean(sampleMhfs)}")
print(f"avg number of participating haps in pop: {np.mean(popHaps)}; avg number of participating haps in samp: {np.mean(sampleHaps)}")
print(f"avg pop marker mut freq in soft sweeps: {np.mean(popMarkerFreqs)}; avg samp marker mut freq in soft sweeps: {np.mean(sampMarkerFreqs)}")
print(f"avg number of gene conversion events per replicate: {np.mean(gcEventCounts)}")
