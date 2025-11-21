#Randomized Motif search
#using Monte Carlo method with repeated random sampling to get a quick, approximate best k-mer to maximize Profile matrix score


import random
from helpers import (
    profile_most_probable_kmer,
    create_profile_matrix,
    score_motifs,
    format_motifs_output
)


def randomized_motif_search(dna_strings, k, t):
    """
    Implement RandomizedMotifSearch algorithm.
    
    Algorithm:
    1. Start with randomly chosen k-mers from each DNA string
    2. Repeat until convergence:
       a. Create profile matrix from current motifs (with pseudocounts)
       b. Find most probable k-mer in each DNA string using the profile
       c. If new motifs have better score, update; otherwise stop
    
    Args:
        dna_strings: List of DNA strings
        k: Length of motifs to find
        t: Number of DNA strings
    
    Returns:
        List of best motifs found in this run
    """

    # Step 1: Randomly select initial motifs
    motifs = []
    for i in range(t):
        # Choose random starting position for k-mer in each DNA string
        max_start = len(dna_strings[i]) - k
        if max_start <= 0:
            # Handle edge case where DNA string is shorter than k
            motifs.append(dna_strings[i][:k].ljust(k, 'A'))
        else:
            start_pos = random.randint(0, max_start)
            motifs.append(dna_strings[i][start_pos:start_pos + k])
    
    best_motifs = motifs[:]
    best_score = score_motifs(best_motifs)
    
    # Step 2: Iteratively improve motifs
    while True:
        # Create profile matrix from current motifs (with pseudocounts)
        profile = create_profile_matrix(motifs)
        
        # Find most probable k-mer in each DNA string
        new_motifs = []
        for i in range(t):
            most_probable = profile_most_probable_kmer(dna_strings[i], k, profile)
            new_motifs.append(most_probable)
        
        # Check if new motifs are better
        new_score = score_motifs(new_motifs)
        
        if new_score < best_score:
            # Better motifs found, continue iteration
            motifs = new_motifs[:]
            best_motifs = new_motifs[:]
            best_score = new_score
        else:
            # No improvement, algorithm has converged
            break
    
    return best_motifs


def run_randomized_motif_search_multiple_times(dna_strings, k, t, num_runs=1000):
    """
    Run RandomizedMotifSearch multiple times and return the best result.
    
    Args:
        dna_strings: List of DNA strings
        k: Length of motifs to find
        t: Number of DNA strings
        num_runs: Number of times to run the algorithm (default: 1000)
    
    Returns:
        Best motifs found across all runs
    """
    best_motifs = None
    best_score = float('inf')
    
    for run in range(num_runs):
        # Run randomized motif search once
        current_motifs = randomized_motif_search(dna_strings, k, t)
        current_score = score_motifs(current_motifs)
        
        # Keep track of best result
        if current_score < best_score:
            best_motifs = current_motifs[:]
            best_score = current_score
        
    
    return best_motifs


def main():
    """
    Read input and run RandomizedMotifSearch 1000 times.
    """
    # Set random seed for reproducibility (optional)
    # random.seed(42)
    
    # Read input file
    with open('data/part_4_data/sample_1.txt', 'r') as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]
    
    # Parse first line: k and t
    k, t = map(int, lines[0].split())
    
    # Parse second line: DNA strings
    dna_strings = lines[1].split()
    
    # Verify we have the expected number of DNA strings
    if len(dna_strings) != t:
        print(f"Error: Expected {t} DNA strings, got {len(dna_strings)}")
        return
    
    # Run RandomizedMotifSearch 1000 times
    print("Running RandomizedMotifSearch 1000 times...", file=__import__('sys').stderr)
    best_motifs = run_randomized_motif_search_multiple_times(dna_strings, k, t, 1000)
    
    # Output result
    result = format_motifs_output(best_motifs)
    print(result)


if __name__ == "__main__":
    main()