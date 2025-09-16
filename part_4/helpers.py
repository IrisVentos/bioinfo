def profile_probability(kmer, profile):
    """
    Calculate the probability of a k-mer given a profile matrix.
    
    Profile matrix format:
    Row 0: A probabilities for each position
    Row 1: C probabilities for each position  
    Row 2: G probabilities for each position
    Row 3: T probabilities for each position
    
    Probability = product of individual position probabilities
    """
    # Map nucleotides to row indices
    nucleotide_to_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    
    probability = 1.0
    
    # Multiply probabilities for each position
    for i, nucleotide in enumerate(kmer):
        if nucleotide in nucleotide_to_index:
            row = nucleotide_to_index[nucleotide]
            probability *= profile[row][i]
        else:
            # Invalid nucleotide - probability becomes 0
            return 0.0
    
    return probability

def profile_most_probable_kmer(text, k, profile):
    """
    Find the k-mer in text with the highest probability according to the profile.
    If there are ties, return the first occurrence.
    """
    max_probability = -1.0
    most_probable_kmer = ""
    
    # Check every possible k-mer in the text
    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]
        probability = profile_probability(kmer, profile)
        
        # Keep track of the highest probability k-mer
        # Note: we only update if probability is STRICTLY greater
        # This ensures we get the first occurrence in case of ties
        if probability > max_probability:
            max_probability = probability
            most_probable_kmer = kmer
    
    return most_probable_kmer


def create_profile_matrix(motifs):
    """
    Create a profile matrix from a collection of k-mers (motifs).
    Uses pseudocounts (Laplace smoothing) by adding 1 to all counts.
    
    Args:
        motifs: List of k-mer strings of equal length
    
    Returns:
        Profile matrix as list of lists [A_row, C_row, G_row, T_row]
        Each row contains probabilities for that nucleotide at each position
    """
    if not motifs or not motifs[0]:
        return []
    
    k = len(motifs[0])  # length of each motif
    t = len(motifs)     # number of motifs
    
    # Initialize count matrix with pseudocounts (1 for each nucleotide at each position)
    # This prevents zero probabilities
    counts = {
        'A': [1] * k,
        'C': [1] * k, 
        'G': [1] * k,
        'T': [1] * k
    }
    
    # Count nucleotides at each position
    for motif in motifs:
        for i, nucleotide in enumerate(motif):
            if nucleotide in counts:
                counts[nucleotide][i] += 1
    
    # Convert counts to probabilities
    # Total count at each position = t + 4 (t motifs + 4 pseudocounts)
    total_count = t + 4
    
    profile = []
    for nucleotide in ['A', 'C', 'G', 'T']:
        prob_row = [count / total_count for count in counts[nucleotide]]
        profile.append(prob_row)
    
    return profile


def score_motifs(motifs):
    """
    Calculate the score of a collection of motifs.
    Score = sum of Hamming distances from each motif to the consensus string.
    Lower score is better (more conserved motifs).
    
    Args:
        motifs: List of k-mer strings of equal length
    
    Returns:
        Integer score (sum of mismatches to consensus)
    """
    if not motifs or not motifs[0]:
        return float('inf')
    
    k = len(motifs[0])
    score = 0
    
    # For each position in the motifs
    for i in range(k):
        # Count nucleotides at this position
        nucleotide_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        
        for motif in motifs:
            if i < len(motif) and motif[i] in nucleotide_counts:
                nucleotide_counts[motif[i]] += 1
        
        # Find the most frequent nucleotide (consensus)
        max_count = max(nucleotide_counts.values())
        
        # Add mismatches to score
        total_motifs = len(motifs)
        mismatches = total_motifs - max_count
        score += mismatches
    
    return score


def format_motifs_output(motifs):
    """
    Format motifs for output as space-separated string.
    """
    return ' '.join(motifs)