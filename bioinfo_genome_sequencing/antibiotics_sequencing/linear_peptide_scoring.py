# Only contiguous substrings, no wrap-around for linear peptides

import sys
from collections import Counter
from peptide_mass_counting import AMINO_ACID_MASS   
from theoretical_spectrum_peptide import linear_spectrum
 
def score(peptide, spectrum):
    """
    Score a linear peptide against a spectrum.
    Score = number of masses in the theoretical spectrum that match
    masses in the experimental spectrum (accounting for multiplicity).
    """
    theoretical = linear_spectrum(peptide)
    
    # Use counters to handle repeated masses correctly
    theo_count = Counter(theoretical)
    exp_count  = Counter(spectrum)
    
    # For each mass, take the minimum count between theoretical and experimental
    total_score = 0
    for mass, count in theo_count.items():
        total_score += min(count, exp_count.get(mass, 0))
    
    return total_score

if __name__ == "__main__":
    input_file = sys.argv[1] if len(sys.argv) > 1 else "bioinfo_genome_sequencing/datasets/dataset_18.txt"
    lines    = open(input_file).read().splitlines()
    peptide  = lines[0].strip()
    spectrum = list(map(int, lines[1].split()))
    print(score(peptide, spectrum))