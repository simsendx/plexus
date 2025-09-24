"""
Oligonucleotide melting temperature calculations.

Based on the oligotm.c library from Primer3.
Copyright (c) 1996-2007 Whitehead Institute for Biomedical Research
"""

import math
from dataclasses import dataclass

# Kelvin to Celsius conversion factor
T_KELVIN = 273.15

class OligotmError(Exception):
    """Base exception for oligotm errors."""
    pass


class InvalidSequenceError(OligotmError):
    """Raised when sequence is invalid or too short."""
    pass


class InvalidConcentrationError(OligotmError):
    """Raised when concentration values are invalid."""
    pass


@dataclass
class TmResult:
    """Result of Tm calculation."""
    Tm: float
    bound: float = 0.0

# SantaLucia (1998) parameters
# Entropy values (DS) in 0.1 cal/(K·mol)
DS_PARAMS = {
    'AA': 222, 'AC': 224, 'AG': 210, 'AT': 204, 'AN': 224,
    'CA': 227, 'CC': 199, 'CG': 272, 'CT': 210, 'CN': 272,
    'GA': 222, 'GC': 244, 'GG': 199, 'GT': 224, 'GN': 244,
    'TA': 213, 'TC': 222, 'TG': 227, 'TT': 222, 'TN': 227,
    'NA': 168, 'NC': 210, 'NG': 220, 'NT': 215, 'NN': 220
}

# Enthalpy values (DH) in 100 cal/mol
DH_PARAMS = {
    'AA': 79, 'AC': 84, 'AG': 78, 'AT': 72, 'AN': 72,
    'CA': 85, 'CC': 80, 'CG': 106, 'CT': 78, 'CN': 78,
    'GA': 82, 'GC': 98, 'GG': 80, 'GT': 84, 'GN': 80,
    'TA': 72, 'TC': 82, 'TG': 85, 'TT': 79, 'TN': 72,
    'NA': 72, 'NC': 80, 'NG': 78, 'NT': 72, 'NN': 72
}

# Delta G values in cal/mol
DG_PARAMS = {
    'AA': 1000, 'AC': 1440, 'AG': 1280, 'AT': 880, 'AN': 880,
    'CA': 1450, 'CC': 1840, 'CG': 2170, 'CT': 1280, 'CN': 1450,
    'GA': 1300, 'GC': 2240, 'GG': 1840, 'GT': 1440, 'GN': 1300,
    'TA': 580, 'TC': 1300, 'TG': 1450, 'TT': 1000, 'TN': 580,
    'NA': 580, 'NC': 1300, 'NG': 1280, 'NT': 880, 'NN': 580
}

def calc_thermodynamics(seq: str) -> tuple:
    """
    Calculate enthalpy (dh) and entropy (ds) for sequence using SantaLucia (1998) parameters.

    Args:
        seq - Sequence to analyze
    """
    dh = 0
    ds = 0
    
    # Select parameter tables
    h_table = DH_PARAMS
    s_table = DS_PARAMS
        
    # Symmetry correction
    if symmetry(seq):
        ds += 14
        
    # Terminal AT penalty
    if seq[0] in 'AT':
        ds += -41
        dh += -23
    elif seq[0] in 'CG':
        ds += 28
        dh += -1
            
    if seq[-1] in 'AT':
        ds += -41
        dh += -23
    elif seq[-1] in 'CG':
        ds += 28
        dh += -1
    
    # Calculate nearest-neighbor contributions
    for i in range(len(seq) - 1):
        pair = seq[i:i+2]
        dh += h_table.get(pair, h_table['NN'])
        ds += s_table.get(pair, s_table['NN'])
    
    return dh, ds


def oligotm(
    seq: str,
    dna_conc: float = 50.0,
    salt_conc: float = 50.0,
    divalent_conc: float = 1.5,
    dntp_conc: float = 0.6,
    dmso_conc: float = 0.0,
    dmso_fact: float = 0.6,
    formamide_conc: float = 0.0,
    annealing_temp: float = 60.0
    ) -> TmResult:
    """
    Calculate melting temperature using SantaLucia nearest-neighbor method.
    
    Args:
        seq: DNA sequence
        dna_conc: DNA concentration (nM)
        salt_conc: Monovalent salt concentration (mM)
        divalent_conc: Divalent cation concentration (mM)
        dntp_conc: dNTP concentration (mM)
        dmso_conc: DMSO concentration (%)
        dmso_fact: DMSO correction factor (default 0.6)
        formamide_conc: Formamide concentration (mol/L)
        annealing_temp: Annealing temperature for binding calculation
    
    Returns:
        TmResult with Tm and bound percentage
        
    Raises:
        InvalidSequenceError: If sequence is too short or contains invalid characters
        InvalidConcentrationError: If concentration values are invalid
    """
    if len(seq) < 2:
        raise InvalidSequenceError("Sequence must be at least 2 bases long")
    
    # Validate sequence contains only valid bases
    valid_bases = set('ATGCN')
    if not all(base in valid_bases for base in seq.upper()):
        raise InvalidSequenceError("Sequence contains invalid characters. Only A, T, G, C, N allowed")

    seq = seq.upper()

    # Validate and convert divalent to monovalent
    dv_to_mv = divalent_to_monovalent(divalent_conc, dntp_conc)

    # Count GC for formamide correction
    gc_count = sum(1 for base in seq if base in 'GC')
    seq_len = len(seq)
    
    # Calculate thermodynamics
    dh, ds = calc_thermodynamics(seq)
    
    delta_H = dh * -100.0  # Convert to cal/mol
    delta_S = ds * -0.1    # Convert to cal/(K·mol)
    
    K_mM = salt_conc + dv_to_mv

    temp = annealing_temp + T_KELVIN
    delta_S = delta_S + 0.368 * (seq_len - 1) * math.log(K_mM / 1000.0)

    ddG = delta_H - temp * delta_S
    ka = math.exp(-ddG / (1.987 * temp))
        
    if symmetry(seq):
        Tm = delta_H / (delta_S + 1.987 * math.log(dna_conc / 1e9)) - T_KELVIN
        bound = (1 / (1 + math.sqrt(1 / ((dna_conc / 1e9) * ka)))) * 100
    else:
        Tm = delta_H / (delta_S + 1.987 * math.log(dna_conc / 4e9)) - T_KELVIN
        bound = (1 / (1 + math.sqrt(1 / ((dna_conc / 4e9) * ka)))) * 100
            
    # Apply DMSO and formamide corrections
    if dmso_conc > 0.0:
        Tm -= dmso_conc * dmso_fact
        Tm += (0.453 * gc_count / seq_len - 2.88) * formamide_conc
    
    return TmResult(Tm=Tm, bound=bound)


def oligodg(seq: str) -> float:
    """
    Calculate delta G of oligonucleotide disruption using SantaLucia parameters.
    
    Raises:
        InvalidSequenceError: If sequence is too short
    """
    if len(seq) < 2:
        raise InvalidSequenceError("Sequence must be at least 2 bases long")
    
    seq = seq.upper()
    dg = -1960  # Initial dG
    
    # Terminal AT penalty
    if seq[0] in 'AT':
        dg += -50
    if seq[-1] in 'AT':
        dg += -50
    
    # Calculate nearest-neighbor contributions
    for i in range(len(seq) - 1):
        pair = seq[i:i+2]
        dg += DG_PARAMS.get(pair, DG_PARAMS['NN'])
    
    if symmetry(seq):
        dg += -430  # Symmetry correction
    
    return dg / 1000.0


# Return the delta G of the last len bases of oligo if oligo is at least len
# bases long; otherwise return the delta G of oligo.
def end_oligodg(seq: str, length: int) -> float:
    """
    Calculate delta G for last 'length' bases of sequence. Return the delta G of the last len bases of oligo if oligo is at least len
   bases long; otherwise return the delta G of oligo.
    
    """
    seq_len = len(seq)
    if seq_len < length:
        return oligodg(seq)
    return oligodg(seq[-length:])


def symmetry(seq: str) -> bool:
    """Check if sequence is self-complementary/symmetrical."""
    seq_len = len(seq)
    if seq_len % 2 == 1:
        return False
    
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    mid = seq_len // 2
    
    for i in range(mid):
        if seq[i] not in complement or seq[-(i+1)] != complement[seq[i]]:
            return False
    
    return True


def divalent_to_monovalent(divalent: float, dntp: float) -> float:
    """Convert divalent salt concentration to monovalent equivalent."""
    if divalent == 0:
        dntp = 0
    if divalent < 0 or dntp < 0:
        raise InvalidConcentrationError("Divalent and dNTP concentrations must be non-negative")
    if divalent < dntp:
        divalent = dntp
    return 120 * math.sqrt(divalent - dntp)

# Both functions return the melting temperature of the given oligo
# calculated as specified by user, but oligotm _should_ only be used on
# DNA sequences of length <= MAX_PRIMER_LENGTH (which is defined
# elsewhere).  seqtm uses oligotm for sequences of length <=
# MAX_PRIMER_LENGTH, and a different, G+C% based formula for longer
# sequences.  For oligotm(), no error is generated on sequences
# longer than MAX_PRIMER_LENGTH, but the formula becomes less
# accurate as the sequence grows longer. 
def long_seq_tm(seq: str,
                start: int,
                length: int,
                salt_conc: float,
                divalent_conc: float,
                dntp_conc: float,
                dmso_conc: float = 0.0,
                dmso_fact: float = 0.6,
                formamide_conc: float = 0.0) -> TmResult:
    """
    Calculate Tm for long sequences using GC% formula.
    
    Uses: Tm = 81.5 + 16.6*log10([Na+]) + 0.41*(%GC) - 600/length
    
    Raises:
        InvalidSequenceError: If sequence bounds are invalid
        InvalidConcentrationError: If concentration values are invalid
    """
    dv_to_mv = divalent_to_monovalent(divalent_conc, dntp_conc)
    salt_conc = salt_conc + dv_to_mv
    
    if start < 0 or length <= 0 or start + length > len(seq):
        raise InvalidSequenceError(f"Invalid sequence bounds: start={start}, length={length}, seq_len={len(seq)}")
    
    subseq = seq[start:start + length].upper()
    gc_count = sum(1 for base in subseq if base in 'GC')
    
    Tm = 81.5 - dmso_conc * dmso_fact + \
         (0.453 * gc_count / length - 2.88) * formamide_conc + \
         16.6 * math.log10(salt_conc / 1000.0) + \
         41.0 * (gc_count / length) - 600.0 / length
    
    return TmResult(Tm=Tm, bound=0.0)

def seqtm(seq: str,
          dna_conc: float = 50.0,
          salt_conc: float = 50.0,
          divalent_conc: float = 1.5,
          dntp_conc: float = 0.6,
          dmso_conc: float = 0.0,
          dmso_fact: float = 0.6,
          formamide_conc: float = 0.0,
          nn_max_len: int = 60,
          annealing_temp: float = 60.0) -> TmResult:
    """
    Calculate melting temperature for sequences of any length.
    
    Uses SantaLucia nearest-neighbor model for sequences <= nn_max_len,
    GC% formula for longer sequences.

    Args:
        nn_max_len - Max Value for which nearest neighbour model is valid. Default is 60.
    """
    seq_len = len(seq)
    
    if seq_len > nn_max_len:
        print("Running longseq")
        return long_seq_tm(seq, 0, seq_len, salt_conc, divalent_conc, dntp_conc,
                          dmso_conc, dmso_fact, formamide_conc)
    else:
        print("Running oligotm")
        return oligotm(seq, dna_conc, salt_conc, divalent_conc, dntp_conc,
                      dmso_conc, dmso_fact, formamide_conc, annealing_temp)

