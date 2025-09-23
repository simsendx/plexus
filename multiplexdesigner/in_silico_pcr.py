# Run in silico PCR to evaluate primers
# https://github.com/pommevilla/ispcr

from ispcr import get_pcr_products

def ispcr_get_producr(primers, sequences):
    """
    get_pcr_products will then iterate through the sequences in sequence_file and find all products amplified by the forward and reverse primer.
    """

    results = get_pcr_products(
        primer_file = primers,
        sequence_file = sequences
    )

    return(results)