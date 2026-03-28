# To compute the theoretical spectrum of a cyclic/linear peptide
# We need the amino acid masses, which we can get from the monoisotopic mass table

# First for linear

AMINO_ACID_MASS = {
    'G': 57,  'A': 71,  'S': 87,  'P': 97,  'V': 99,
    'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
    'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131,
    'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186
}

def linear_spectrum(peptide):
    n = len(peptide)
    prefix_mass = [0] * (n + 1)
    for i in range(1, n + 1):
        aa = peptide[i - 1]
        prefix_mass[i] = prefix_mass[i - 1] + AMINO_ACID_MASS[aa]

    spectrum = [0]
    for i in range(n):
        for j in range(i + 1, n + 1):
            spectrum.append(prefix_mass[j] - prefix_mass[i])

    return sorted(spectrum)

with open("bioinfo_genome_sequencing/datasets/dataset_14.txt") as f:
    peptide = f.read().strip()

result = linear_spectrum(peptide)
#print(" ".join(map(str, result)))

def cyclic_spectrum(peptide):
    n = len(peptide)
    prefix_mass = [0] * (n + 1)
    for i in range(1, n + 1):
        aa = peptide[i - 1]
        prefix_mass[i] = prefix_mass[i - 1] + AMINO_ACID_MASS[aa]

    peptide_mass = prefix_mass[n]

    spectrum = [0]
    for i in range(n):
        for j in range(i + 1, n + 1):
            fragment = prefix_mass[j] - prefix_mass[i]
            spectrum.append(fragment)
            # Add wraparound subpeptide (cyclic case), skip full peptide and empty
            if i > 0 and j < n:
                spectrum.append(peptide_mass - fragment)

    return sorted(spectrum)

with open("bioinfo_genome_sequencing/datasets/dataset_15.txt") as f:
    peptide = f.read().strip()

result = cyclic_spectrum(peptide)
print(" ".join(map(str, result)))