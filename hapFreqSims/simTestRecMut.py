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
            #print(f"sample genotypes for {n} individuals: {n}")
            genotypes = []
        else:
            genomes = self.genomes
        for genome in genomes:
            num += genome.chromosome1.hasMarkerMut() + genome.chromosome2.hasMarkerMut()
            if n:
                genotypes.append((genome.chromosome1.markerMut, genome.chromosome1.selMut))
                genotypes.append((genome.chromosome2.markerMut, genome.chromosome2.selMut))
            denom += 2
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
        return self.selMut != 0
    
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

def doConversionAtSelSite(donorChrom, recipChrom):
    if donorChrom.selMut > 0 and recipChrom.selMut == 0:
        recipChrom.population.currSelMutCount += 1
        recipChrom.selMut = recipChrom.population.currSelMutCount
        
def makeRecombinantGametesForParent(parent, totalGcRate):
    newParent = parent.clone()
    if random.random() < totalGcRate:
        if random.randint(0, 1):
            doConversionAtSelSite(newParent.chromosome2, newParent.chromosome1)
        else:
            doConversionAtSelSite(newParent.chromosome1, newParent.chromosome2)
    
    return newParent.chromosome1, newParent.chromosome2

def makeRecombinantOffspringChroms(parent1, parent2, totalGcRate):
    parent1Gametes = makeRecombinantGametesForParent(parent1, totalGcRate)
    parent1Gamete = random.choice(parent1Gametes)
    
    parent2Gametes = makeRecombinantGametesForParent(parent2, totalGcRate)
    parent2Gamete = random.choice(parent2Gametes)
    
    return parent1Gamete, parent2Gamete
        
def makeOffspringGenome(parent1, parent2, totalGcRate, recMutRate):
    offspringChrom1, offspringChrom2 = makeRecombinantOffspringChroms(parent1, parent2, totalGcRate)
    
    # possibility of recurrent mutation to the adaptive allele
    if random.random() < recMutRate:
        offspringChrom1.addMutation()
    if random.random() < recMutRate:
        offspringChrom2.addMutation()
        
    offspringGenome = Genome(offspringChrom1, offspringChrom2)
    return offspringGenome

def runSimUntilFixation(Ngen, Nt, pop, totalGcRate, recMutRate):
    done = 0
    popH21, sampleH21 = None, None
    established = False
    while not done:
        for gen in range(Ngen):
            currN = Nt[gen]
            genomes = []

            for offspring in range(currN):
                parent1 = pop.selectParent()
                parent2 = pop.selectParent()
                offspringGenome = makeOffspringGenome(parent1, parent2, totalGcRate, recMutRate)
                genomes.append(offspringGenome)
            pop = Population(genomes=genomes, s=pop.s, h=pop.h)
            pop.updateFitnesses()
            selMutFreq = pop.getSelMutFreq()
            if selMutFreq == 0:
                #print(f"selected mutation lost at generation: {gen}. restarting")
                #init again:
                pop = Population(N=Nanc, s=pop.s, h=pop.h)
                targetGenomeIndex = pop.selectRandomGenomeIndex()
                pop.genomes[targetGenomeIndex].addSelMutAndRemoveMarker()
                pop.updateFitnessOfIndexedGenome(targetGenomeIndex)
                break
            elif (not established) and selMutFreq > 0.1:
                if verbose:
                    print(f"selected mutation established at generation: {gen}")
                established=True
            if selMutFreq == 1:
                if verbose:
                    print(f"selected mutation fixed at generation: {gen}")
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
    return gen, popH21, popMhf, sampleH21, sampleMhf, popNumParticipatingHaps, sampleNumParticipatingHaps

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

baseNanc, gcScalar, recMutScalar, s, h, Q, numReps = sys.argv[1:]
verbose=False

gcScalar = float(gcScalar)
recMutScalar = float(recMutScalar)
s=float(s)
h=float(h)
Q=float(Q)
numReps = int(numReps)

baseNanc=float(baseNanc)
Nanc = int(baseNanc/Q)
Ngen = Nanc*10
totalGcRate = gcScalar*5e-06*Q
recMutRate = recMutScalar*1e-08*Q
sys.stderr.write(f"running {numReps} reps with N={baseNanc}, s={s}, h={h}, 4Ng={4*Nanc*totalGcRate} and 4Nu={4*Nanc*recMutRate}\n")
s = s*Q


gens, popH21s, popMhfs, sampleH21s, sampleMhfs, popHaps, sampleHaps = [], [], [], [], [], [], []

for rep in range(numReps):
    print(f"starting rep {rep}")
    gen, popH21, popMhf, sampleH21, sampleMhf, popNumParticipatingHaps, sampleNumParticipatingHaps = runReplicate(
        Nanc, Ngen, totalGcRate, recMutRate, s, h)
    print(f"done with rep {rep}")

    gens.append(gen)
    if popH21 != 0:
        popH21s.append(popH21)
        popMhfs.append(popMhf)
    if sampleH21 != 0:
        sampleH21s.append(sampleH21)
        sampleMhfs.append(sampleMhf)
    popHaps.append(popNumParticipatingHaps)
    sampleHaps.append(sampleNumParticipatingHaps)

    print("")

popSofts = [x for x in popHaps if x > 1]
sampleSofts = [x for x in sampleHaps if x > 1]

print(f"fraction of soft sweeps in pop: {len(popSofts)/len(popHaps)}; fraction of soft sweeps in sample: {len(sampleSofts)/len(sampleHaps)}")
print(f"avg pop h2/h1 of soft sweeps: {np.mean(popH21s)}; avg sample h2/h1 of soft sweeps: {np.mean(sampleH21s)}")
print(f"avg minor hap pop freq of soft sweeps: {np.mean(popMhfs)}; avg minor hap sample freq of soft sweeps: {np.mean(sampleMhfs)}")
print(f"avg number of participating haps in pop: {np.mean(popHaps)}; avg number of participating haps in samp: {np.mean(sampleHaps)}")
