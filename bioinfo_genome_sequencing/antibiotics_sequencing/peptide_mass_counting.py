# Using dynamic programming to count the number of linear peptides of a given mass m, where the mass is determined by the amino acid composition. 
# We will use the integer masses of the 20 standard amino acids, but since I and L have the same mass (113) and K and Q have the same mass (128), we will only consider 18 unique masses.


AMINO_ACID_MASS = {
    'G': 57,  'A': 71,  'S': 87,  'P': 97,  'V': 99,
    'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
    'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131,
    'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186
}

# Use only the 18 unique masses (I/L and K/Q collapsed)
UNIQUE_MASSES = list(set(AMINO_ACID_MASS.values()))  # 18 distinct values

def count_linear_peptides(m):
    """
    Count linear peptides of integer mass m using dynamic programming.
    dp[x] = number of linear peptides with total mass x
    """
    dp = [0] * (m + 1)
    dp[0] = 1  # base case: one empty peptide has mass 0

    for mass in range(1, m + 1):
        for aa_mass in UNIQUE_MASSES:
            if mass - aa_mass >= 0:
                dp[mass] += dp[mass - aa_mass]

    return dp[m]

m = 1429
print(count_linear_peptides(m))
