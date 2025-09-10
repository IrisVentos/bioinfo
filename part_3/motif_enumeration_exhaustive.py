#could reuse function neighbors in part_2 with dictionary of neighborhood
#let's rewind from beginning

import sys
import os
sys.path.append('part_2')
from helpers import hamming_distance, neighbors

def get_all_kmers(dna_string, k):
    """Get all k-mers from a DNA string."""
    kmers = []
    for i in range(len(dna_string) - k + 1):
        kmers.append(dna_string[i:i + k])
    return kmers

def appears_in_dna_with_mismatches(pattern, dna_string, d):
    """Check if pattern appears in dna_string with at most d mismatches."""
    k = len(pattern)
    for i in range(len(dna_string) - k + 1):
        kmer = dna_string[i:i + k]
        if hamming_distance(pattern, kmer) <= d:
            return True
    return False

def motif_enumeration(dna_list, k, d):
    """Find all (k,d)-motifs in DNA strings."""
    patterns = set()
    
    # For each k-mer in each DNA string
    for dna_string in dna_list:
        kmers = get_all_kmers(dna_string, k)
        for kmer in kmers:
            # Generate all neighbors within d mismatches
            pattern_neighbors = neighbors(kmer, d)
            
            # Check each neighbor
            for pattern in pattern_neighbors:
                # Check if this pattern appears in ALL DNA strings with at most d mismatches
                appears_in_all = True
                for dna in dna_list:
                    if not appears_in_dna_with_mismatches(pattern, dna, d):
                        appears_in_all = False
                        break
                
                if appears_in_all:
                    patterns.add(pattern)
    
    return sorted(list(patterns))

def main():
    # Read input file
    input_file = "data/part_3_data/sample_1.txt"
    
    try:
        with open(input_file, 'r') as f:
            lines = f.read().strip().split('\n')
        
        # Parse first line: k and d
        k, d = map(int, lines[0].split())
        
        # Parse second line: DNA strings
        dna_strings = lines[1].split()
        
        print(f"k = {k}, d = {d}")
        print(f"DNA strings: {dna_strings}")
        print()
        
        # Run motif enumeration
        result = motif_enumeration(dna_strings, k, d)
        
        # Output results
        print("Found motifs:")
        print(' '.join(result))
        
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        print("Please create the file with the following format:")
        print("Line 1: k d")
        print("Line 2: DNA strings separated by spaces")
        print()
        print("Example content for motif_input.txt:")
        print("3 1")
        print("ATTTGGC TGCCTTA CGGTATC GAAAATT")

if __name__ == "__main__":
    main()