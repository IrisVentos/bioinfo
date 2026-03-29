from theoretical_spectrum_peptide import AMINO_ACID_MASS, linear_spectrum, cyclic_spectrum

UNIQUE_MASSES = sorted(set(AMINO_ACID_MASS.values()))

def is_consistent(peptide, spectrum):
    peptide_spectrum = linear_spectrum(peptide)
    spectrum_copy = list(spectrum)
    for m in peptide_spectrum:
        if m in spectrum_copy:
            spectrum_copy.remove(m)
        else:
            return False
    return True

def expand(candidates):
    return [peptide + [aa] for peptide in candidates for aa in UNIQUE_MASSES]

def all_rotations(peptide):
    return [tuple(peptide[i:] + peptide[:i]) for i in range(len(peptide))]

def cyclopeptide_sequencing(spectrum):
    spectrum = sorted(spectrum)
    target_mass = max(spectrum)

    candidates = [[]]
    final_peptides = []
    seen = set()

    while candidates:
        candidates = expand(candidates)
        next_candidates = []

        for peptide in candidates:
            m = sum(peptide)

            if m == target_mass:
                if cyclic_spectrum(peptide) == spectrum:
                    # check all rotations to avoid duplicates
                    rotations = all_rotations(peptide)
                    if not any(r in seen for r in rotations):
                        for r in rotations:
                            seen.add(r)
                        final_peptides.append(peptide)
            elif is_consistent(peptide, spectrum):
                next_candidates.append(peptide)

        candidates = next_candidates

    return final_peptides

with open("bioinfo_genome_sequencing/datasets/dataset_16.txt") as f:
    spectrum = list(map(int, f.read().strip().split()))

results = cyclopeptide_sequencing(spectrum)

# each peptide dash-separated, peptides space-separated
print(" ".join("-".join(map(str, p)) for p in results))