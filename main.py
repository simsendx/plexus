import os
import primer3
#import pysam
from Bio import SeqIO
from src.hello import say_hello_to

# Get mutation list


# Get metadata
min_primer_length = 18
max_primer_lgenth = 30
min_amplicon_size = 40
max_amplicon_size = 80
min_padding_bases = 3


# Import design regions

# For each mutation
#   get sequence +/- 200 bp
#   define design regions
    

# Use pysam to get region from fasta
# https://pysam.readthedocs.io/en/latest/api.html


# ---------------// Hairpin analysis //------------------
# Use primer3 to calculate hairpin
# max seq length is 60 nucleotides
# uses the thal library from primer3
# results stored as ThermoResult object
hairpin_thermo_result = primer3.bindings.calc_hairpin(
    seq = "ACTCTGTCTGTGTG",  # DNA sequence to analyze for hairpins
    mv_conc = 50,            # Monovalent cation conc. (mM)
    dv_conc = 1.5,           # Divalent cation conc. (mM)
    dntp_conc = 0.6,         # total dNTP conc. (mM)
    dna_conc = 50.0,         # DNA conc. (nM)
    temp_c = 37.0,           # Simulation temperature for dG (Celsius)
    max_loop = 30,           # Maximum loop size in structures
    output_structure = False
)

print(hairpin_thermo_result)

say_hello_to("Stefan")