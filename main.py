import os
import primer3
#import pysam
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from src.hello import say_hello_to
import pandas as pd

say_hello_to("Stefan")

def calculate_primer_design_regions(j_min, j_max, design_parameters):
    """
    Calculate the forward and reverse primer design regions based on junction positions
    and design parameters.
    
    Args:
        j_min (int): Left position of the junction
        j_max (int): Right position of the junction
        design_parameters (dict): Dictionary containing design parameters including:
            - min_primer_length (int): Minimum primer length
            - max_primer_length (int): Maximum primer length
            - min_amplicon_gap (int): Minimum amplicon gap
            - max_amplicon_gap (int): Maximum amplicon gap
            - padding_bases (int): Number of padding bases
            
    Returns:
        tuple: Two lists containing [start, end] positions for left and right primer regions
    """
    
    # Extract parameters
    max_primer_length = design_parameters['max_primer_length']
    max_amplicon_gap = design_parameters['max_amplicon_gap']
    padding_bases = design_parameters['padding_bases']
    
    # Calculate forward primer region
    fp_start = j_max - max_primer_length - max_amplicon_gap + padding_bases
    fp_end = j_min + padding_bases
    
    # Calculate reverse primer region
    rp_start = j_max + padding_bases
    rp_end = j_min + max_primer_length + max_amplicon_gap - padding_bases
    
    # Return both regions as a tuple
    return [fp_start, fp_end], [rp_start, rp_end]

def depth_first_search_with_pruning(node, depth, alpha, beta, maximizing_player):
    # Base case: reached maximum depth or leaf node
    if depth == 0 or is_leaf_node(node):
        return evaluate_node(node)
    
    if maximizing_player:
        max_eval = float('-inf')
        for child in get_children(node):
            eval = depth_first_search_with_pruning(child, depth - 1, alpha, beta, False)
            max_eval = max(max_eval, eval)
            alpha = max(alpha, eval)
            
            # Pruning step: beta cut-off
            if beta <= alpha:
                break  # Cut off this branch
        return max_eval
    else:
        min_eval = float('inf')
        for child in get_children(node):
            eval = depth_first_search_with_pruning(child, depth - 1, alpha, beta, True)
            min_eval = min(min_eval, eval)
            beta = min(beta, eval)
            
            # Pruning step: alpha cut-off
            if beta <= alpha:
                break  # Cut off this branch
        return min_eval
    
# TODO how to handle overlapping junctions? There can
# be different cases:
# (1) Junctions very close together -> merge immediately
# (2) Moderately far apart junctions -> allow up to 14 nt 5' overlap
# (3) Identify junctions which cannot be merged -> report and/or ignore one of them

def read_junction_coordinates(csv_path):
    """
    Read junction coordinates from CSV file and convert to list of pairs.
    
    Args:
        csv_path (str): Path to the CSV file containing junction coordinates
        
    Returns:
        list: List of [min, max] coordinate pairs
    """
    df = pd.read_csv(csv_path)
    
    junctions = []
    for _, row in df.iterrows():
        five_prime = int(row['Five_Prime_Coordinate'])
        three_prime = int(row['Three_Prime_Coordinate'])
        junctions.append([min(five_prime, three_prime), max(five_prime, three_prime)])
    
    return junctions

def sort_junctions(junctions):
    """
    Sort junctions by their start position.
    
    Args:
        junctions (list): List of [min, max] coordinate pairs
        
    Returns:
        list: Sorted list of junction pairs
    """
    return sorted(junctions, key=lambda x: x[0])

def merge_adjacent_junctions(sorted_junctions, max_gap):
    """
    Merge junctions that are within the specified maximum gap.
    
    Args:
        sorted_junctions (list): Sorted list of [min, max] coordinate pairs
        max_gap (int): Maximum allowed gap between junctions to merge
        
    Returns:
        list: List of merged junction coordinates
    """
    if not sorted_junctions:
        return []
        
    merged_junctions = []
    current = sorted_junctions[0]
    
    for next_junction in sorted_junctions[1:]:
        if current[1] + max_gap >= next_junction[0]:
            # Merge the junctions
            current = [
                min(current[0], next_junction[0]),
                max(current[1], next_junction[1])
            ]
        else:
            merged_junctions.append(current)
            current = next_junction
    
    merged_junctions.append(current)
    return merged_junctions

def process_and_merge_junctions(csv_path, design_parameters):
    """
    Import junctions from CSV file, process them, and merge where appropriate.
    
    Args:
        csv_path (str): Path to the CSV file containing junction coordinates
        design_parameters (dict): Dictionary containing design parameters including max_amplicon_gap
        
    Returns:
        list: List of merged junction coordinates as [min, max] pairs
    """
    # Read and convert coordinates
    junctions = read_junction_coordinates(csv_path)
    
    # Sort junctions
    sorted_junctions = sort_junctions(junctions)
    
    # Merge adjacent junctions
    merged_junctions = merge_adjacent_junctions(
        sorted_junctions, 
        design_parameters['max_amplicon_gap']
    )
    
    return merged_junctions

# ---------------// Hairpin analysis //------------------
# Use primer3 to calculate hairpin
# max seq length is 60 nucleotides
# uses the thal library from primer3
# results stored as ThermoResult object

def calculate_hairpin(sequence, pcr_conditions):
    """
    Calculate hairpin thermodynamics for a DNA sequence.
    
    Args:
        sequence (str or Seq): DNA sequence to analyze
        pcr_conditions (dict): PCR conditions including ion concentrations
        
    Returns:
        ThermoResult: Primer3 thermodynamic analysis result
    """
    return primer3.bindings.calc_hairpin(
        seq=str(sequence),
        mv_conc=pcr_conditions['mv_concentration'],
        dv_conc=pcr_conditions['dv_concentration'],
        dntp_conc=pcr_conditions['dntp_concentration'],
        dna_conc=pcr_conditions['dna_concentration'],
        temp_c=37.0,  # Standard temperature for dG calculations
        max_loop=30,
        output_structure=False
    )

def calculate_melting_temperature(sequence, pcr_conditions):
    """
    Calculate melting temperature for a DNA sequence.
    
    Args:
        sequence (str or Seq): DNA sequence to analyze
        pcr_conditions (dict): PCR conditions including concentrations and temperature
        
    Returns:
        float: Melting temperature in Celsius
    """
    return primer3.bindings.calc_tm(
        seq=str(sequence),
        mv_conc=pcr_conditions['mv_concentration'],
        dv_conc=pcr_conditions['dv_concentration'],
        dntp_conc=pcr_conditions['dntp_concentration'],
        dna_conc=pcr_conditions['dna_concentration'],
        dmso_conc=pcr_conditions['dmso_concentration'],
        dmso_fact=pcr_conditions['dmso_fact'],
        formamide_conc=pcr_conditions['formamide_concentration'],
        annealing_temp_c=pcr_conditions['annealing_temperature'],
        max_nn_length=60,
        tm_method='santalucia',
        salt_corrections_method='santalucia'
    )

def calculate_homodimer(sequence, pcr_conditions):
    """
    Calculate homodimer formation thermodynamics.
    
    Args:
        sequence (str or Seq): DNA sequence to analyze
        pcr_conditions (dict): PCR conditions including ion concentrations
        
    Returns:
        ThermoResult: Primer3 thermodynamic analysis result
    """
    return primer3.bindings.calc_homodimer(
        str(sequence),
        mv_conc=pcr_conditions['mv_concentration'],
        dv_conc=pcr_conditions['dv_concentration'],
        dntp_conc=pcr_conditions['dntp_concentration'],
        dna_conc=pcr_conditions['dna_concentration'],
        temp_c=37.0,
        max_loop=30,
        output_structure=False
    )

def calculate_heterodimer(sequence1, sequence2, pcr_conditions):
    """
    Calculate heterodimer formation thermodynamics between two sequences.
    
    Args:
        sequence1 (str or Seq): First DNA sequence
        sequence2 (str or Seq): Second DNA sequence
        pcr_conditions (dict): PCR conditions including ion concentrations
        
    Returns:
        ThermoResult: Primer3 thermodynamic analysis result
    """
    return primer3.bindings.calc_heterodimer(
        str(sequence1),
        str(sequence2),
        mv_conc=pcr_conditions['mv_concentration'],
        dv_conc=pcr_conditions['dv_concentration'],
        dntp_conc=pcr_conditions['dntp_concentration'],
        dna_conc=pcr_conditions['dna_concentration'],
        temp_c=37.0,
        max_loop=30,
        output_structure=False
    )

# Example usage:
# candidate_sequence = Seq("TTGTGAGTTTTTGAAATCTCTGTGA")
# hairpin_result = calculate_hairpin(candidate_sequence, pcr_conditions)
# tm = calculate_melting_temperature(candidate_sequence, pcr_conditions)
# homodimer_result = calculate_homodimer(candidate_sequence, pcr_conditions)
# heterodimer_result = calculate_heterodimer(
#     candidate_sequence, 
#     candidate_sequence.reverse_complement(), 
#     pcr_conditions
# )

#------// Workflow //------
# 1) Prepare inputs
# 2) Design single-plex solutions for each target
# 3) Check the specificity of each target
# 4) Pick the optimal multiplex solution
# 5) Output multiplex

#-----// (1) Prepare inputs //---------
# Get the mutation list (junctions) around which to design primer pairs

# Target genome
panel_name  = "test_panel"
genome      = "hg38"
sample_type = "cfDNA"
protocol    = "simsen"

# Get metadata
design_parameters = {
    "primer_length_range": [18, 30],
    "amplicon_gap_range": [20, 80],
    "junction_padding_bases": 3,
    "initial_solutions": 100,
    "top_solutions_to_keep": 4,
    "minimum_plexity": 5,
    "maximum_plexity": 20,
    "target_plexity": 20,
    "force_plexity": False,  # If True, the target_plexity will be used to force the plexity of the primers
    "variant_threshold": 0.01
}

pcr_conditions = {
    "annealing_temperature": 60,
    "primer_concentration": 50,
    "dntp_concentration": 0.6,
    "dna_concentration": 50,
    "mv_concentration": 50,
    "dv_concentration": 1.5,
    "dmso_concentration": 0.0,
    "dmso_fact": 0.6,
    "formamide_concentration": 0.8
}

penalties = {
    "snp_penalty": 1.0,
    "primer_length_penalty": 1.0,
    "primer_complexity_penalty": 1.0,
    "polyA_penalty": 5,
    "polyT_penalty": 5,
    "polyC_penalty": 10,
    "polyG_penalty": 10,
    "amplicon_length_penalty": 1.0
}

# Define and merge junctions
# Each junction has a le
# ft (j_min) and right (j_max) position
j_min = 123456
j_max = 123456

left_primer_region, right_primer_region = calculate_primer_design_regions(j_min, j_max, design_parameters)


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
