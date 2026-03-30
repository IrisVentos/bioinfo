# Adapting cycloptide sequencing for spectra with noise and missing peaks
# In general, if a mass occurs m times in the theoretical spectrum of Peptide and n times in the experimental spectrum Spectrum, 
# then it contributes the minimum of m and n to Score(Peptide, Spectrum)

# Score(Peptide, Spectrum) = sum over all masses of min(m, n) = columns shared by the two spectra in the spectrum graph


import sys
from collections import Counter
from peptide_mass_counting import AMINO_ACID_MASS   
from theoretical_spectrum_peptide import cyclic_spectrum
 
def score(peptide, spectrum):
    """
    Score a cyclic peptide against a spectrum.
    Score = number of masses in the theoretical spectrum that match
    masses in the experimental spectrum (accounting for multiplicity).
    """
    theoretical = cyclic_spectrum(peptide)
    
    # Use counters to handle repeated masses correctly
    theo_count = Counter(theoretical)
    exp_count  = Counter(spectrum)
    
    # For each mass, take the minimum count between theoretical and experimental
    total_score = 0
    for mass, count in theo_count.items():
        total_score += min(count, exp_count.get(mass, 0))
    
    return total_score

if __name__ == "__main__":
    input_file = sys.argv[1] if len(sys.argv) > 1 else "bioinfo_genome_sequencing/datasets/dataset_17.txt"
    lines    = open(input_file).read().splitlines()
    peptide  = lines[0].strip()
    spectrum = list(map(int, lines[1].split()))
    print(score(peptide, spectrum))