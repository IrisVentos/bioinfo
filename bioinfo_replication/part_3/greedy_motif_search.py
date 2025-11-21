from helpers import (
    profile_most_probable_kmer,
    create_profile_matrix, 
    score_motifs,
    format_motifs_output
)


def greedy_motif_search(dna_strings, k, t):
    """
    Implement the GreedyMotifSearch algorithm.
    
    Algorithm:
    1. Initialize BestMotifs with first k-mer from each string
    2. For each possible k-mer in the first string:
       a. Set this as Motif1
       b. Build profile from Motif1
       c. Find profile-most probable k-mer in string 2 → Motif2
       d. Build profile from Motif1, Motif2
       e. Find profile-most probable k-mer in string 3 → Motif3
       f. Continue until all strings processed
       g. Score the resulting motif collection
       h. Update BestMotifs if this collection scores better
    3. Return BestMotifs
    
    Args:
        dna_strings: List of DNA sequence strings
        k: Length of motifs to find
        t: Number of strings (should equal len(dna_strings))
    
    Returns:
        List of k-mers representing the best motif collection
    """
    # Step 1: Initialize BestMotifs with first k-mer from each string
    best_motifs = []
    for dna_string in dna_strings:
        if len(dna_string) >= k:
            best_motifs.append(dna_string[:k])
        else:
            # Handle edge case: string shorter than k
            best_motifs.append(dna_string)
    
    best_score = score_motifs(best_motifs)
    
    print(f"Initial BestMotifs: {best_motifs}")
    print(f"Initial score: {best_score}")
    print()
    
    # Step 2: Try each k-mer in the first string as Motif1
    first_string = dna_strings[0]
    
    for i in range(len(first_string) - k + 1):
        # Step 2a: Set current k-mer as Motif1
        motif1 = first_string[i:i + k]
        current_motifs = [motif1]
        
        print(f"Trying Motif1 = '{motif1}' (position {i})")
        
        # Step 2b-e: Greedily find best k-mers for remaining strings
        for j in range(1, t):  # j goes from 1 to t-1 (strings 2 to t)
            # Build profile from current motifs
            profile = create_profile_matrix(current_motifs)
            
            # Find profile-most probable k-mer in string j+1
            best_kmer = profile_most_probable_kmer(dna_strings[j], k, profile)
            current_motifs.append(best_kmer)
            
            print(f"  Added Motif{j+1} = '{best_kmer}' from string {j+1}")
        
        # Step 2f-g: Score the current motif collection
        current_score = score_motifs(current_motifs)
        print(f"  Final motifs: {current_motifs}")
        print(f"  Score: {current_score}")
        
        # Step 2h: Update BestMotifs if current collection is better
        if current_score < best_score:
            best_motifs = current_motifs[:]  # Make a copy
            best_score = current_score
            print(f"  ✓ New best motifs found!")
        else:
            print(f"  No improvement (best score still {best_score})")
        
        print()
    
    print(f"=== FINAL RESULT ===")
    print(f"Best motifs: {best_motifs}")
    print(f"Best score: {best_score}")
    
    return best_motifs


def display_motif_alignment(motifs):
    """
    Display the motif alignment nicely formatted.
    """
    print("\nMotif Alignment:")
    for i, motif in enumerate(motifs):
        print(f"String {i+1}: {motif}")
    
    # Show consensus
    if motifs and motifs[0]:
        k = len(motifs[0])
        consensus = ""
        for pos in range(k):
            # Count nucleotides at this position
            counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
            for motif in motifs:
                if pos < len(motif) and motif[pos] in counts:
                    counts[motif[pos]] += 1
            
            # Find most common nucleotide
            most_common = max(counts.items(), key=lambda x: x[1])[0]
            consensus += most_common
        
        print(f"Consensus: {consensus}")


def main():
    # Read input from file
    try:
        with open('data/part_3_data/sample_4.txt', 'r') as file:
            lines = file.read().strip().split('\n')
    except FileNotFoundError:
        print("Error: input.txt file not found!")
        print("\nExpected format:")
        print("Line 1: k t (space-separated integers)")
        print("Line 2: space-separated DNA strings")
        return
    
    if len(lines) < 2:
        print("Error: Input file must have at least 2 lines")
        return
    
    # Parse input
    try:
        k, t = map(int, lines[0].split())
        dna_strings = lines[1].split()
    except ValueError:
        print("Error: Could not parse input. Check format.")
        return
    
    if len(dna_strings) != t:
        print(f"Error: Expected {t} DNA strings, but got {len(dna_strings)}")
        return
    
    print(f"GreedyMotifSearch Parameters:")
    print(f"k = {k} (motif length)")
    print(f"t = {t} (number of strings)")
    print(f"DNA strings: {dna_strings}")
    print()
    
    # Run GreedyMotifSearch
    best_motifs = greedy_motif_search(dna_strings, k, t)
    
    # Display results
    display_motif_alignment(best_motifs)
    
    # Output in required format
    print(f"\n=== OUTPUT ===")
    print(format_motifs_output(best_motifs))


if __name__ == "__main__":
    main()