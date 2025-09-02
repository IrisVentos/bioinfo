# The goal is to find all k-mers with up to d mismatches with Pattern in Text 

from helpers import hamming_distance
from itertools import product

def neighbors(pattern, d):
    """
    Generate all k-mers within Hamming distance d from pattern.
    This is the d-neighborhood of pattern.
    
    Args:
        pattern (str): The original k-mer pattern
        d (int): Maximum Hamming distance allowed
        
    Returns:
        set: All k-mers within distance d from pattern
    """
    if d == 0:
        return {pattern}
    
    if len(pattern) == 1:
        return {'A', 'C', 'G', 'T'}
    
    neighborhood = set()
    suffix_neighbors = neighbors(pattern[1:], d)
    
    for suffix in suffix_neighbors:
        if hamming_distance(pattern[1:], suffix) < d:
            # If suffix uses < d mismatches, first nucleotide can be anything
            for nucleotide in ['A', 'C', 'G', 'T']:
                neighborhood.add(nucleotide + suffix)
        else:
            # If suffix uses exactly d mismatches, first nucleotide must match
            neighborhood.add(pattern[0] + suffix)
    
    return neighborhood

def frequent_words_with_mismatches(text, k, d):
    """
    Find all most frequent k-mers with up to d mismatches in text.
    
    Args:
        text (str): The DNA string to analyze
        k (int): Length of k-mers to consider
        d (int): Maximum number of mismatches allowed
        
    Returns:
        list: All most frequent k-mers with up to d mismatches
    """
    patterns = []
    freq_map = {}  # Dictionary to count approximate matches
    n = len(text)
    
    # For each k-mer in the text
    for i in range(n - k + 1):
        pattern = text[i:i + k]  # Extract k-mer at position i
        
        # Get all k-mers within distance d of this pattern
        neighborhood = neighbors(pattern, d)
        
        # Increment count for each neighbor
        for neighbor in neighborhood:
            if neighbor in freq_map:
                freq_map[neighbor] += 1
            else:
                freq_map[neighbor] = 1
    
    # Find the maximum frequency
    max_count = max(freq_map.values()) if freq_map else 0
    
    # Collect all k-mers with maximum frequency
    for pattern in freq_map:
        if freq_map[pattern] == max_count:
            patterns.append(pattern)
    
    return patterns

# Read input from file
with open('data/part_2_data/sample_6.txt', 'r') as f:
    lines = f.read().strip().split('\n')
    text = lines[0]
    k = int(lines[1])
    d = int(lines[2])

# Find frequent words with mismatches
result = frequent_words_with_mismatches(text, k, d)

# Print results
print(' '.join(result))
