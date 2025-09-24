import primer3
import math
from itertools import combinations

class ThermoObject:
    """Class to represent primer3 output"""
    
    def __init__(self, junction_name: str, chrom: str, five_prime: int, three_prime: int):
        self.junction_name = junction_name
        self.chrom = chrom
        self.five_prime = five_prime
        self.three_prime = three_prime
        self.design_region = None
        self.design_start = None
        self.design_end = None
        self.gc_content = None
        self.tm = None
        self.hairpin = None
        self.homodimer = None


    def primer3_calculate_oligo(
        self, 
        primer_list, 
        min_gc: float = 30.0, 
        max_gc: float = 70.0, 
        min_tm: float = 57.0, 
        max_tm: float = 63.0,
        primer_max_self_any_th: float = 45.0,
        primer_max_hairpin_th: float = 24.0
        ):

        oligo_calc = primer3.thermoanalysis.ThermoAnalysis()

        oligo_calc.set_thermo_args(
            mv_conc = 50.0,
            dv_conc = 1.5,
            dntp_conc = 0.6,
            dna_conc = 50.0,
            dmso_conc = 0.0,
            dmso_fact = 0.6,
            formamide_conc = 0.0,
            annealing_temp_c = -10,
            temp_c = 37.0,
            max_nn_length = 60,
            max_loop = 30,
            tm_method = 'santalucia',
            salt_corrections_method = 'santalucia'
        )

        for primer in primer_list:

            # Check GC content of primer
            primer_gc = gc_content(primer) 
            if primer_gc < min_gc | primer_gc > max_gc:
                continue
            
            # Check primer melting temperature
            primer_tm = oligo_calc.calc_tm(primer)
            if primer_tm < min_tm | primer_tm > max_tm:
                continue

            # Check thermodynamic alignment of oligo
            primer_homodimer = oligo_calc.calc_homodimer(primer)
            if primer_homodimer.tm > primer_max_self_any_th:
                continue

            primer_hairpin = oligo_calc.calc_hairpin(primer)
            if primer_hairpin.tm > primer_max_hairpin_th:
                continue

        for a, b in combinations(primer_list, 2):
            print(f"Primer 1: {a} + Primer 2: {b}")
            print(oligo_calc.calc_heterodimer(a, b))
            print(oligo_calc.calc_end_stability(a, b))


    
    # ---------------// Hairpin analysis //------------------
    # Use primer3 to calculate hairpin
    # max seq length is 60 nucleotides
    # uses the thal library from primer3
    # results stored as ThermoResult object

    def primer3_calculate_hairpin(sequence, pcr_conditions):
        """
        Calculate hairpin thermodynamics for a DNA sequence.
        
        Args:
            sequence (str or Seq): DNA sequence to analyze
            pcr_conditions (dict): PCR conditions including ion concentrations
            
        Returns:
            ThermoResult: Primer3 thermodynamic analysis result
        """
        return primer3.bindings.calc_hairpin(
            seq = str(sequence),
            mv_conc = pcr_conditions['mv_concentration'],
            dv_conc = pcr_conditions['dv_concentration'],
            dntp_conc = pcr_conditions['dntp_concentration'],
            dna_conc = pcr_conditions['dna_concentration'],
            temp_c = 37.0,  # Standard temperature for dG calculations
            max_loop = 30,
            output_structure = False
        )


def calculate_amount_bound(delta_H, delta_S, K_mM, DNA_nM, divalent_conc, dntp_conc, length, annealing_temp, sym, salt_correction = "santalucia1998", T_KELVIN = 273.15):

    # Based on primer3:
    # https://github.com/libnano/primer3-py/blob/2242f67811b89aeeccca878b55767d954204c3c5/primer3/src/libprimer3/oligotm.c#L740
    if salt_correction == "santalucia1998":
        K_mM = K_mM + divalent_to_monovalent(divalent_conc, dntp_conc)
        delta_S = delta_S + 0.368 * (length - 1) * math.log(K_mM / 1000.0)
        if sym == 1:
            ddG = delta_H - (annealing_temp + T_KELVIN) * delta_S
            ka = math.exp(-ddG / (1.987 * (annealing_temp + T_KELVIN)))
            percent_bound = (1 / (1 + math.sqrt(1 / ((DNA_nM / 1e9) * ka)))) * 100
        else:
            ddG = delta_H - (annealing_temp + T_KELVIN) * delta_S
            ka = math.exp(-ddG / (1.987 * (annealing_temp + T_KELVIN)))
            percent_bound = (1 / (1 + math.sqrt(1/((DNA_nM/4000000000.0) * ka)))) * 100

    return(percent_bound)

def divalent_to_monovalent(divalent: float, dntp: float) -> float:
    """
    Convert divalent salt concentration to monovalent equivalent.
    
    Parameters:
        divalent (float): concentration of divalent cations
        dntp (float): concentration of dNTPs (nucleotides)
    
    Returns:
        float: equivalent monovalent concentration
    """
    if divalent == 0:
        dntp = 0
    if divalent < 0 or dntp < 0:
        raise ValueError("divalent and dNTP concentrations must be non-negative")
    if divalent < dntp:
        # Melting temp does not depend on divalent cations
        divalent = dntp
    return 120.0 * math.sqrt(divalent - dntp)


def symmetry(seq: str) -> int:
    """
    Return 1 if DNA sequence is symmetrical, 0 otherwise.
    Symmetry means base pairing: A<->T, C<->G.
    """
    seq_len = len(seq)
    
    # Odd-length sequences cannot be symmetrical
    if seq_len % 2 == 1:
        return 0

    # Define complement rules
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # Compare each character from start and end
    for i in range(seq_len // 2):
        if seq[i] not in complement or complement[seq[i]] != seq[-(i + 1)]:
            return 0

    return 1

def gc_content(sequence: str) -> float:
    """
    Calculate the GC content of a DNA sequence.
    
    Parameters:
        sequence (str): DNA sequence consisting of A, T, G, C.
    
    Returns:
        float: GC content as a percentage.
    """
    if not sequence:
        return 0.0  # Handle empty string safely
    
    sequence = sequence.upper()  # Ensure case-insensitivity
    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100

def primer3_calculate_melting_temperature(sequence, pcr_conditions):
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

def primer3_calculate_homodimer(sequence, pcr_conditions):
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

def primer3_calculate_heterodimer(sequence1, sequence2, pcr_conditions):
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

def primer3_calculate_end_stability():
    """
    A ThermoResult object with thermodynamic characteristics of the 3’ hybridization interaction
    """
    return primer3.bindings.calc_end_stability()


def primer3_design_primers(seq_args, global_args, misprime_lib = None, mishyb_lib = None):
    """
    Run the Primer3 design process.

    Args:
        seq_args (Dict[str, Any]): Primer3 sequence/design args as per Primer3 docs
        global_args (Dict[str, Any]): Primer3 global args as per Primer3 docs
        misprime_lib (Optional[Dict[str, Any]]): Sequence name: sequence dictionary for mispriming checks.
        mishyb_lib (Optional[Dict[str, Any]]): Sequence name: sequence dictionary for mishybridization checks.

    Returns:
        A dictionary of Primer3 results
    """
    return primer3.bindings.design_primers(seq_args, global_args, misprime_lib, mishyb_lib)