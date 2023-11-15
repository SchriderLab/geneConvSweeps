# This monstrosity was copied and heavily edited from some code that Logan
# Whitehouse wrote for the Timesweeper method (github.com/SchriderLab/Timesweeper)

import sys, os, re
import numpy as np


def get_slim_code(slim_file):
    with open(slim_file, "r") as infile:
        raw_lines = [line.replace("\n", "") for line in infile.readlines()]

    return raw_lines

def get_slim_info(slim_lines):
    def get_ints(lines, query):
        """Splits a line by regex split and grabs all digits. Returns as list."""
        # Uses max to handle max_years_b0, but might need to check if that causes other parsing issues
        return [
            int(s)
            for s in re.split("\W+", [i for i in lines if query in i][0])
            if s.isdigit()
        ]

    def get_float(lines, query):
        for line in lines:
            if query in line:
                return float(line.split('"generation_time", ')[1].split(")")[0])

    # format: off
    Q = get_ints(slim_lines, """defineConstant("Q",""")[0]
    gen_time = get_float(slim_lines, """defineConstant("generation_time",""")
    burn_time_mult = get_ints(slim_lines, """defineConstant("burn_in",""")[0]
    max_years_b0 = max(get_ints(slim_lines, """defineConstant("_T","""))
    physLen = int(slim_lines[slim_lines.index("    _recombination_ends = c(")+1].strip().split(")")[0])+1
    # format: on

    pop_sizes_line = [i for i in slim_lines if "_N" in i][0]
    pop_sizes_start = slim_lines.index(pop_sizes_line) + 2
    # The 4,3 will need to be more flexible
    pop_sizes_end = (
        slim_lines[pop_sizes_start:].index(
            '    defineConstant("num_epochs", length(_T));'
        )
        - 2
        + pop_sizes_start
    )  # End of definition idx, will be end of sampling ep array

    sizes = []
    for i in slim_lines[pop_sizes_start:pop_sizes_end]:
        sizes.append([int(s) for s in re.split("\W+", i) if s.isdigit()])
    first_biggest = max([i[0] for i in sizes])
    burn_in_gens = first_biggest * burn_time_mult

    return Q, gen_time, max_years_b0, round(burn_in_gens), physLen


def inject_sweep_type(raw_lines, sweep, selCoeff, domCoeff, selTime, popNumber, chromLen, dumpFilePrefix):
    if sweep in ["hard"]:
        raw_lines[raw_lines.index("    defineConstant(\"burn_in\", 0.0);")] = "    defineGlobal(\"burn_in\", 0.0);"
        raw_lines.insert(
            raw_lines.index("    initializeMutationRate(Q*c("),
            f"    defineConstant('sweep', '{sweep}');\n    defineConstant('dump_file_name', \"{dumpFilePrefix}.\" + rep_number + \".dump\");",
        )
        raw_lines.insert(
            raw_lines.index("    initializeMutationRate(Q*c("),
            f"    defineConstant('selCoeff', Q * {selCoeff});",
        )
        raw_lines.insert(
            raw_lines.index("    // `5.5 Rescaling population sizes to improve simulation performance`.")+1,
            f"    defineConstant('sweepSite', asInteger(round({chromLen}/2)));",
        )
        raw_lines.insert(
            raw_lines.index("    initializeMutationRate(Q*c("),
            f"    defineConstant('targetPop', {popNumber});\n",
        )
        raw_lines.insert(
            raw_lines.index("    initializeMutationRate(Q*c("),
            f"    defineGlobal('fixedBeforeEnd', F);\n",
        )
        raw_lines.insert(
            raw_lines.index("    community.registerLateEvent(NULL, \"{dbg(self.source); end();}\", G_end, G_end);"),
            f"    community.registerLateEvent(NULL, \"{{checkOnSweep();}}\", selTime, G_end);\n    community.registerLateEvent(NULL, \"{{endSimBlock();}}\", G_end, G_end);\n",
        )
        raw_lines.insert(
            raw_lines.index("    G = G_start + gdiff(T_start, _T);")+1,
            f"""    Ge = max(G);
    ts={selTime};
    N0 = max(N[0, targetPop], N[0, 0]);
    sTime = Ge - asInteger(N0*4*ts);

    if (sTime < 0)
    {{
        totalBurninGens = -sTime + 1;
        defineGlobal("burn_in", totalBurninGens / N_max);
        G_start = community.tick + asInteger(round(burn_in * N_max));
        T_start = max(_T);
        G = G_start + gdiff(T_start, _T);
        Ge = max(G);
        sTime = Ge - asInteger(N0*4*ts);
    }}
    defineConstant("G_end", Ge);
    cat('Total sim time: ' + G_end + ' generations\\n');

    if (sTime < 1) {{cat("Warning: sTime < 1 after adding burn-in: " + sTime + "\\n"); sTime = 1;}}

    defineConstant('selTime', sTime);

    cat('Sweep starts at generation ' + selTime + '\\n');
""")
        raw_lines.pop(raw_lines.index("    G_end = max(G);"))
        raw_lines.insert( #insert our selected mutation (m2) and our marker to find recombinants (m3)
            raw_lines.index("    initializeMutationRate(Q*c("),
            f"    initializeMutationType('m2', {domCoeff}, 'f', selCoeff);\n    initializeMutationType('m3', 0.0, 'f', 0.0);\n    m3.convertToSubstitution = F;",
        )

def inject_gene_conv(raw_lines, gc_ratio):
    ratesLineIndex = raw_lines.index("    _recombination_rates = c(")
    raw_lines[ratesLineIndex] = f"    _recombination_rates = 500*{gc_ratio}*c(0,"
    raw_lines[ratesLineIndex+1] = raw_lines[ratesLineIndex+1].replace(");", ",0);")
    raw_lines[raw_lines.index("    _recombination_ends = c(")] = f"    _recombination_ends = c(sweepSite-1,sweepSite,"
    raw_lines.insert(
            raw_lines.index("    defineConstant(\"recombination_rates\", _recombination_rates);")+1,
            f"    initializeGeneConversion(1.0, 2, 1.0);"
        )

def sanitize_slim(raw_lines):
    new_lines = []
    sanitized_lines = 0
    mode = 1

    for line in raw_lines:
        unwantedLine = ("sim.treeSeqOutput(trees_file, metadata=metadata);" in line and mode == 2) or \
                   "initializeTreeSeq();" in line or \
                   """"inds=p"+pop+".sampleIndividuals("+n+"); " +""" in line or \
                   "defineConstant(\"trees_file\"" in line
        if not unwantedLine:
            new_lines.append(line)
        else:
            sanitized_lines += 1

        if mode == 1:
            if line == "function (void)end(void) {":
                mode = 2
        if mode == 2:
            if "sim.treeSeqOutput(trees_file, metadata=metadata);" in line:
                mode = 3

    return new_lines, sanitized_lines


def inject_sampling(raw_lines, pop, samp_counts, gens, outfile_path):
    samp_eps_line = [i for i in raw_lines if "sampling_episodes" in i][0]
    samp_eps_start = raw_lines.index(
        samp_eps_line
    )  # Start of sampling_episodes constant idx
    samp_eps_end = (
        raw_lines[samp_eps_start:].index("    ), c(3, 3)));") + samp_eps_start
    )  # End of definition idx, will be end of sampling ep array

    new_lines = []
    new_lines.extend(raw_lines[samp_eps_start : samp_eps_start + 1])
    last_line_id = len(samp_counts)
    for samps, gen, idx in zip(samp_counts, gens, range(last_line_id)):
        if idx < last_line_id - 1:
            new_lines.append("        " + f"c({pop[1:]}, {samps}, {gen}),")
        else:
            new_lines.append("        " + f"c({pop[1:]}, {samps}, {gen})")
    new_lines.append("    " + f"), c(3, {len(gens)})));")

    for line in raw_lines:
        if "treeSeqRememberIndividuals" in line:
            raw_lines[
                raw_lines.index(line)
            ] = f"""            "{pop}.outputSample("+n+", replace=F, filePath='{outfile_path}', append=T);}}","""

    finished_lines = raw_lines[:samp_eps_start]
    finished_lines.extend(new_lines)
    finished_lines.extend(raw_lines[samp_eps_end + 1 :])

    return finished_lines


def make_sel_blocks(pop, sample_size):
    intro_blocks = f"""\n1 late() {{
    sim.subpopulations.individuals.genomes.addNewDrawnMutation(m3, 1);
}}

s1 1 late(){{
    // save the state of the simulation
    sim.outputFull(dump_file_name);
    
    // introduce the sweep mutation
    // Note: this code below may not work if we have more than 2 subpops!
    if (length(sim.subpopulations) > targetPop)
    {{
        defineConstant("mutIntroPop", targetPop);
    }}
    else
    {{
        defineConstant("mutIntroPop", 0);
    }}
    defineGlobal("currTargetPop", mutIntroPop);
    target = sample(sim.subpopulations[currTargetPop].genomes, 1);
    target.removeMutations();
    target.addNewDrawnMutation(m2, sweepSite);
}}

"""

    check_blocks = f"""function (void)endSimBlock(void) {{
        if (!fixedBeforeEnd)
        {{
            cat("NO FIXATION BEFORE SIMULATION END!\\n");
            sim.subpopulations[currTargetPop].outputMSSample(200, replace=F, filterMonomorphic=F);
            d = deleteFile(dump_file_name);
            cat("Simulation finished. Attempted to delete dump file. Success?: " + d + "\\n");
            sim.simulationFinished();
        }}
}}

function (void)checkOnSweep(void) {{
    fixed = 0;
    polymorphic = 1;

    if (sim.countOfMutationsOfType(m2) == 0)
    {{
        if (sum(sim.substitutions.mutationType == m2) == 1)
        {{
            fixed = 1;
        }}
        polymorphic = 0;
    }}
    else
    {{
        if (sim.mutationFrequencies(sim.subpopulations[currTargetPop], sim.mutationsOfType(m2))[0] == 1.0)
        {{
            fixed = 1;
            polymorphic = 0;
        }}
    }}
    
    if (fixed)
    {{
        cat("FIXED at generation " + community.tick + "\\n");
        sim.subpopulations[currTargetPop].outputMSSample(200, replace=F, filterMonomorphic=F);
        m3freq = sum(sim.mutationFrequencies(sim.subpopulations[currTargetPop], sim.mutationsOfType(m3)));
        cat("total m3 frequency: " + m3freq + "\\n");
        d = deleteFile(dump_file_name);
        cat("Simulation finished. Attempted to delete dump file. Success?: " + d + "\\n");
        defineGlobal('fixedBeforeEnd', T);
        sim.simulationFinished();
    }}
    else if (!polymorphic)
    {{
        cat("LOST - RESTARTING\\n");

        // go back to when the sweep was introduced
        sim.readFromPopulationFile(dump_file_name);

        // start a newly seeded run
        setSeed(rdunif(1, 0, asInteger(2^62) - 1));

        // re-introduce the sweep mutation
        defineGlobal("currTargetPop", mutIntroPop);
        target = sample(sim.subpopulations[currTargetPop].genomes, 1);
        target.removeMutations();
        target.addNewDrawnMutation(m2, sweepSite);
    }}
}}

s2 1:2 late() {{
    if(currTargetPop != targetPop)
    {{
        if(length(sim.subpopulations) > targetPop)
        {{
            defineGlobal("currTargetPop", targetPop);
        }}
    }}
}}

"""


    reschedule_block = f"""1 early(){{
    community.rescheduleScriptBlock(s1, start=selTime, end=selTime);
    community.rescheduleScriptBlock(s2, start=selTime, end=G_end);
}}

"""

    all_blocks = []
    for block in intro_blocks, check_blocks, reschedule_block:
        all_blocks.extend(block.split("\n"))

    return all_blocks

def removeOrigSamplingBlock(raw_lines):
    new_lines = []
    skipMode = 0
    skipCount = 0
    for line in raw_lines:
        if skipMode == 0 and "    // Sample individuals." in line:
            skipMode = 1
        if skipMode == 1:
            skipCount += 1
            if line == "    }":
                skipMode = 2
        else:
            new_lines.append(line)
    return new_lines

def write_slim(finished_lines, out_file_name):
    with open(out_file_name, "w") as outfile:
        for line in finished_lines:
            outfile.write(line + "\n")

def main():
    scriptFileName, dumpFilePrefix, outFileName = sys.argv[1:]
    params = scriptFileName.split("/")[-1].rstrip(".slim")
    specName, demogModel, mapName, contigName, targetPop, popNumber, sampleSize, Q, selCoeff, domCoeff, selTime, geneConvScalar = params.split("-")
    popNumber, sampleSize = [int(x) for x in [popNumber, sampleSize]]
    Q, selCoeff, domCoeff, selTime, geneConvScalar = [float(x) for x in [Q, selCoeff, domCoeff, selTime, geneConvScalar]]

    # Info scraping and calculations
    raw_lines = get_slim_code(scriptFileName)
    raw_lines[raw_lines.index(" * stdpopsim 0.2.0")] = " * stdpopsim 0.2.0 -- with modifications by Logan Whitehouse and Dan Schrider's weird pipeline"
    raw_lines, sanitized_line_count = sanitize_slim(raw_lines)

    Q, gen_time, max_years_b0, burn_in_gens, physLen = get_slim_info(raw_lines)
    burn_in_gens = int(round(burn_in_gens / Q))
    burn_in_years = burn_in_gens * gen_time

    print(max_years_b0, Q, burn_in_years, gen_time)
    end_gen = int(
        round(((max_years_b0 / Q) + burn_in_years) / gen_time)
    )  # Convert from earliest year from bp to gens, note that max_years_b0 is not scaled by Q

    # Logging - is there a cleaner way to do this?
    print("geneConvSweeps SLiM Injection")
    print("Q Scaling Value:", Q)
    print("Gen Time:", gen_time)
    print("Simulated Chrom Length:", physLen)
    print()
    print("Burn in years:", burn_in_gens * gen_time)
    print("Burn in gens:", burn_in_gens)
    print()
    print("Max number gens simulated (post-burn):", end_gen)
    print()
    print(
        "Number Years Simulated (inc. burn):", max_years_b0 + (burn_in_gens * gen_time)
    )
    print("Number gens simulated (inc. burn):", end_gen)
    print()
    print("Selection start in 4N gen ago:", selTime)
    print()
    print(f"Sampled pop number: p{popNumber}")
    print()
    print("Sample size:", sampleSize)
    print()
    print(f"Removed {sanitized_line_count} lines in sanitize_slim")

    # Injection
    inject_sweep_type(raw_lines, "hard", selCoeff, domCoeff, selTime, popNumber, physLen, dumpFilePrefix)

    inject_gene_conv(raw_lines, geneConvScalar)

    selection_lines = make_sel_blocks(popNumber, sampleSize)

    raw_lines = removeOrigSamplingBlock(raw_lines)
    finished_lines = []
    finished_lines.extend(raw_lines)
    finished_lines.extend(selection_lines)

    write_slim(finished_lines, outFileName)

    print("Done!")
    print("Output written to:", outFileName)


if __name__ == "__main__":
    main()

