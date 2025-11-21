#big difference with greedy_motif_search function is that it takes into consideration the unfair scoring and inaccurate oversimplification of P(rare event)=0
#Cromwell's rule + Laplace's rule of succession (substitutions of zeroes with small numbers called pseudocounts)

"""
GreedyMotifSearch with Pseudocounts Implementation

Code Challenge: Implement GreedyMotifSearch with pseudocounts.
Input: Integers k and t, followed by a space-separated collection of strings Dna.
Output: A collection of strings BestMotifs resulting from applying 
        GreedyMotifSearch(Dna, k, t) with pseudocounts (+1 in count motifs matrix).
"""

from helpers import (
    profile_most_probable_kmer,
    create_profile_matrix,
    score_motifs,
    format_motifs_output
)


def greedy_motif_search(dna, k, t):
    """
    Greedy algorithm to find motifs with pseudocounts.
    
    Algorithm:
    1. Initialize BestMotifs as first k-mer from each string
    2. For each k-mer in the first string:
        a. Set it as first motif
        b. For each remaining string, find most probable k-mer given current profile
        c. If resulting motifs have better score, update BestMotifs
    
    Args:
        dna: List of DNA strings
        k: Length of motifs to find
        t: Number of strings (should equal len(dna))
    
    Returns:
        List of best motifs found
    """
    # Initialize BestMotifs with first k-mer from each string
    best_motifs = [string[:k] for string in dna]
    best_score = score_motifs(best_motifs)
    
    # Try each k-mer from the first string as starting point
    for i in range(len(dna[0]) - k + 1):
        # Start with i-th k-mer from first string
        motifs = [dna[0][i:i + k]]
        
        # For each remaining string, find the most probable k-mer
        for j in range(1, t):
            # Create profile from current motifs (with pseudocounts)
            profile = create_profile_matrix(motifs)
            
            # Find most probable k-mer in string j given the profile
            most_probable = profile_most_probable_kmer(dna[j], k, profile)
            motifs.append(most_probable)
        
        # Check if these motifs are better than current best
        current_score = score_motifs(motifs)
        if current_score < best_score:
            best_motifs = motifs[:]
            best_score = current_score
    
    return best_motifs


def main():
    """
    Main function to read input and run GreedyMotifSearch.
    """
    # Read input from file or stdin
    try:
        with open('data/part_3_data/sample_5.txt', 'r') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
    except FileNotFoundError:
        print("Input data missing")
        print("Format:")
        print("Line 1: k t (space-separated integers)")
        print("Line 2: DNA strings (space-separated)")
        return
    
    if len(lines) < 2:
        print("Error: Input file must contain at least 2 lines")
        return
    
    # Parse first line: k and t
    try:
        k, t = map(int, lines[0].split())
    except ValueError:
        print("Error: First line must contain two integers k and t")
        return
    
    # Parse second line: DNA strings
    dna_strings = lines[1].split()
    
    if len(dna_strings) != t:
        print(f"Error: Expected {t} DNA strings, got {len(dna_strings)}")
        return
    
    # Validate that all strings are long enough
    for i, dna_string in enumerate(dna_strings):
        if len(dna_string) < k:
            print(f"Error: DNA string {i+1} is shorter than k={k}")
            return
    
    # Run GreedyMotifSearch
    best_motifs = greedy_motif_search(dna_strings, k, t)
    
    # Output results
    result = format_motifs_output(best_motifs)
    print(result)
    


if __name__ == "__main__":
    main()