from collections import defaultdict

def create_profile(motifs):
    """Create profile matrix from motifs with pseudocounts"""
    k = len(motifs[0])  # length of motifs
    profile = defaultdict(lambda: [1] * k)  # pseudocount of 1
    
    # Count occurrences
    for motif in motifs:
        for i, nucleotide in enumerate(motif):
            profile[nucleotide][i] += 1
    
    # Convert to probabilities
    total = len(motifs) + 4  # number of motifs + 4 pseudocounts
    for nucleotide in ['A', 'C', 'G', 'T']:
        for i in range(k):
            profile[nucleotide][i] /= total
    
    return profile

def get_probability(kmer, profile):
    """Calculate probability of kmer given profile"""
    prob = 1.0
    for i, nucleotide in enumerate(kmer):
        prob *= profile[nucleotide][i]
    return prob

def get_all_kmers(dna, k):
    """Get all k-mers from a DNA string"""
    kmers = []
    for i in range(len(dna) - k + 1):
        kmers.append(dna[i:i+k])
    return kmers

def motifs_from_profile(profile, dna_list):
    """Find most probable k-mer in each DNA string given profile"""
    k = len(profile['A'])
    motifs = []
    
    for dna in dna_list:
        best_kmer = ""
        best_prob = 0
        
        for kmer in get_all_kmers(dna, k):
            prob = get_probability(kmer, profile)
            if prob > best_prob:
                best_prob = prob
                best_kmer = kmer
        
        motifs.append(best_kmer)
    
    return motifs

# Given data
dna = [
    "TGACGTTC",
    "TAAGAGTT", 
    "GGACGAAA",
    "CTGTTCGC"
]

initial_motifs = ["TGA", "GTT", "GAA", "TGT"]

print("Initial motifs:", initial_motifs)
print()

# Step 1: Create profile from initial motifs
profile = create_profile(initial_motifs)

print("Profile matrix:")
for nucleotide in ['A', 'C', 'G', 'T']:
    print(f"{nucleotide}: {[round(p, 3) for p in profile[nucleotide]]}")
print()

# Step 2: Find most probable 3-mer in each DNA string
new_motifs = motifs_from_profile(profile, dna)

print("After one iteration:")
print("New motifs:", new_motifs)
print()
print("Answer:", " ".join(new_motifs))