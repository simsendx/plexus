import os
import primer3
#import pysam
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from src.hello import say_hello_to

say_hello_to("Stefan")

#------// Workflow //------
# Prepare inputs

# Design single-plex solutions for each target

# Check the specificity of each target

# Pick the optimal multiplex solution

# Output multiplex

#-----// (1) Prepare inputs //---------
# Get mutation list

# Target genome
panel_name = "test_panel"
genome = "hg38"
sample_type = "cfDNA"
protocol = "simsen"

# Get metadata
min_primer_length = 18
max_primer_length = 30
min_amplicon_gap = 25
max_amplicon_gap = 80
padding_bases = 3

# Define and merge junctions
j_min = 123456
j_max = 123456

# Define primer design regions based on junctions
fp_start = j_max - max_primer_length - max_amplicon_gap + padding_bases
fp_end = j_min + padding_bases

rp_start = j_max + padding_bases
rp_end = j_min + max_primer_length + max_amplicon_gap - padding_bases

left_primer_region = [fp_start, fp_end]
right_primer_region = [rp_start, rp_end]

# TODO how to handle overlapping junctions? There can
# be different cases:
# (1) Junctions very close together -> merge immediately
# (2) Moderately far aprt junctions -> allow up to 14 nt 5' overlap
# (3) Identify junctions which cannot be merged -> report and/or ignore one of them

# Use pysam to get region from fasta
# https://pysam.readthedocs.io/en/latest/api.html


# Define "tail-sequences" or adapters to be appended to each FP or RP, respectively

fp_adapter = "forward adapter"
rp_adapter = "reverse adapter"

"""
>hg38_dna range=chr3:179195014-179195414 5'pad=200 3'pad=200 strand=+ repeatMasking=none
AAAAAAAAAAAAAAAAACACAGGGTAGAGACAGATAACAAGGGATTCCTG
ACACCTACGCTGTAGATAGCTATAACATTTCAATAGGAATCTTGGGAATC
TGGGGAACTGGAAGCAGCTGGGAAAGCACCGGCATCTGAGCCACTCAGTA
CTCTTCCAGGCTTGTGTGTGAGCCTTGCCACCTTCAGGTATTAGCACTTG
AAATCTAACTTCTTTATGAAGCTCCTTATTTACTTGCCTTCTCGGTGAAA
AAAAAAAACAACAACAAAAACCGTTTCCTTCCCCATCTGGTGCCAGCACT
GTATAATTCCTAATAAGCTTGAAAAGATAATGTTGGCAATACTTTACAAG
TTTATTGCTAATGTTTAACATTTATTAGGTGCTTGCTGTGTGCTGATCTC
T
"""


#-----// (2) Designer //---------

# DESIGNER is the core engine for creating all the primer and probe candidates and computing all the
# thermodynamic aspects of design such as target unimolecular folding, primer folding, bimolecular
# hybridization, and solving the multi-state coupled equilibria for the amount bound for the desired bimolecular
# duplex. The amount bound is directly related to the amount of signal generated in a diagnostic assay.
# DESIGNER also analyzes each primer and probe design for a series of heuristic properties such as
# sequence complexity, polyG test, oligo length penalty, amplicon length penalty, etc. Each of the scoring effects
# are multiplied by weighting factors and combined into an overall score for each primer/ probe set.

# In December 2018 functionality was added to PanelPlex is support for dbSNP (for human genome
# version GRCh38 and mouse genome version mm10). PanelPlex automatically detects all the positions with
# alternative alleles with CAF > 0.01 within each of the design regions and then automatically designs the
# primers to AVOID those SNP sites.

# Hybridization conditions


# Designer parameters

init_solutions = 100
top_solutions_to_keep = 4

# 1. Generate potential solutions TODO Generate k-mers as initial solutions
# 2. Score each solutions TODO define scoring algorithm, including SNP penalty
# 3. Select top candidates for each target (init_solutions)
# 4. BLAST against amplicons in panel and discard matches TODO implement BLAST search
# 5. Keep top N solutions (top_solutions_to_keep)
# 6. Perform multiplex picking in N^X space for N solutions and X targets TODO implement multiplex picker algorithm


#-----// Generate k-mer solutions for targets

region = "CACAGGGTAGAGACAGATAACAAGGGATTCCTGACACCAAAAAAAAAAATACGCTGTAGATAGCTATAACATTTCAATAGGAATCTTGGGAATC"

k_min = min_primer_length
k_max = max_primer_length

def generate_kmers(k_min, k_max, sequence):
    kmers = []
    for k in range (k_min, k_max):
        # max position to search
        max_pos = len(sequence) + 1 - k
        
        for x in range(max_pos):
            kmer = sequence[x:x+k]
            kmers.append(kmer)

    return(kmers)


obj = generate_kmers(k_min, k_max, region)
print(obj)
print(len(obj))

# ---------------// Hairpin analysis //------------------
# Use primer3 to calculate hairpin
# max seq length is 60 nucleotides
# uses the thal library from primer3
# results stored as ThermoResult object


candidate_sequence = Seq("TTGTGAGTTTTTGAAATCTCTGTGA")

hairpin_thermo_result = primer3.bindings.calc_hairpin(
    seq = str(candidate_sequence),  # DNA sequence to analyze for hairpins
    mv_conc = 50,            # Monovalent cation conc. (mM)
    dv_conc = 1.5,           # Divalent cation conc. (mM)
    dntp_conc = 0.6,         # total dNTP conc. (mM)
    dna_conc = 50.0,         # DNA conc. (nM)
    temp_c = 37.0,           # Simulation temperature for dG (Celsius)
    max_loop = 30,           # Maximum loop size in structures
    output_structure = False
)

print(hairpin_thermo_result)

melting_temp = primer3.bindings.calc_tm(
    seq = str(candidate_sequence), 
    mv_conc = 50.0, 
    dv_conc = 1.5, 
    dntp_conc = 0.6, 
    dna_conc = 50.0, 
    dmso_conc = 0.0,
    dmso_fact = 0.6,
    formamide_conc = 0.8,
    annealing_temp_c = 60,
    max_nn_length = 60,
    tm_method='santalucia',
    salt_corrections_method='santalucia'
)

print(melting_temp)

homo_dimer = primer3.bindings.calc_homodimer(
    str(candidate_sequence), 
    mv_conc = 50.0, 
    dv_conc = 1.5, 
    dntp_conc = 0.6, 
    dna_conc = 50.0, 
    temp_c = 37.0, 
    max_loop = 30, 
    output_structure = False
)

print(homo_dimer)

het_dimer = primer3.bindings.calc_heterodimer(
    str(candidate_sequence), 
    str(candidate_sequence.reverse_complement()), 
    mv_conc = 50.0, 
    dv_conc = 1.5, 
    dntp_conc = 0.6, 
    dna_conc = 50.0, 
    temp_c = 37.0, 
    max_loop = 30, 
    output_structure = False
)

print(het_dimer)