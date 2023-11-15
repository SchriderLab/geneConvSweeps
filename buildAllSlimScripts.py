import os
import stdpopsim
from contextlib import redirect_stdout

speciesToQ = {"HomSap": 10, "DroMel": 100, "AraTha": 10}
speciesLs = ["HomSap", "DroMel", "AraTha"]
sampleSizes = [200]
defaultModel = "const1pop"
demogModels = {"HomSap": ["OutOfAfrica_2T12", defaultModel], "DroMel": ["African3Epoch_1S16", defaultModel], "AraTha": ["African3Epoch_1H18", defaultModel]}
targetPops = {("HomSap", "OutOfAfrica_2T12"): ["AFR", "EUR"], ("HomSap", defaultModel): [0], \
              ("DroMel", "African3Epoch_1S16"): ["AFR"], ("DroMel", defaultModel): [0], \
              ("AraTha", "African3Epoch_1H18"): ["SouthMiddleAtlas"], ("AraTha", defaultModel): [0]}
speciesToContig = {"HomSap": "chr12", "DroMel": "chr3L", "AraTha": "chr4"}
defaultMap = "flat"
speciesToMap = {"HomSap": [defaultMap], "DroMel": [defaultMap], "AraTha": [defaultMap]}
selCoeffs = [0.1, 0.01, 0.001]
domCoeffs = [0.0, 0.5, 1.0]
selTimes = [0.01, 0.05, 0.1, 0.5, 1]
geneConvRatios = [0, 1, 5, 10, 50]
scriptOutDir = "slimScripts"
os.system(f"mkdir -p {scriptOutDir}")

engine = stdpopsim.get_engine('slim')
for specName in speciesLs:
    species = stdpopsim.get_species(specName)
    contigName = speciesToContig[specName]

    for demogModel in demogModels[specName]:
        if demogModel == defaultModel:
            t = int(round(4*species.population_size)) #have to trick it to run for ~4N generations without a burn-in
            model = stdpopsim.PiecewiseConstantSize(species.population_size, (t, species.population_size))
        else:
            model = species.get_demographic_model(demogModel)
        [pop.id for pop in model.populations]

        for targetPop in targetPops[(specName, demogModel)]:
            for sampleSize in sampleSizes:
                if demogModel == defaultModel:
                    samplesRequested = [sampleSize]
                    assert targetPop == 0
                    popNumber = 0
                else:
                    samplesRequested = []
                    popNumber = None
                    for i in range(len(model.populations)):
                        pop = model.populations[i]
                        if pop.name == targetPop:
                            samplesRequested.append(sampleSize)
                            popNumber = i
                        else:
                            samplesRequested.append(0)
                    assert popNumber != None
                print(demogModel, demogModel==defaultModel, targetPop, popNumber)

                #not actually using this! we will inject our own sampling lines after the sweeping mut fixes
                #samples = model.get_samples(*samplesRequested)
                samples = model.get_samples(*[0]*len(samplesRequested))

                for mapName in speciesToMap[specName]:
                    if mapName != defaultMap:
                        genetic_map = mapName
                    else:
                        genetic_map = None
                contig = species.get_contig(contigName, genetic_map=genetic_map, length_multiplier=0.001)

                Q = speciesToQ[specName]
                for selCoeff, domCoeff, selTime in [(x1, x2, x3) for x1 in selCoeffs for x2 in domCoeffs for x3 in selTimes]:
                    for geneConvRatio in geneConvRatios:
                        scriptFileName = scriptOutDir + "/" + f"{specName}-{demogModel}-{mapName}-{contigName}-{targetPop}-{popNumber}-{sampleSize}-{Q}-{selCoeff}-{domCoeff}-{selTime}-{geneConvRatio}.slim"
                        with open(scriptFileName, 'w') as f:
                            with redirect_stdout(f):
                                ts = engine.simulate(model, contig, samples, slim_script=True, verbosity=2, slim_scaling_factor=Q, slim_burn_in=0)
