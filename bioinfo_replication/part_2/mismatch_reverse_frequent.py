# most frequent_words_with_mismatches_and_reverse_complements

from helpers import hamming_distance
from d_neighborhood import neighbors

def reverse_complement(pattern):
    """
    Generate the reverse complement of a DNA string.
    
    Args:
        pattern (str): DNA string
        
    Returns:
        str: Reverse complement of the pattern
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement[nucleotide] for nucleotide in pattern[::-1])

def approximate_pattern_count(text, pattern, d):
    """
    Count approximate occurrences of pattern in text with up to d mismatches.
    
    Args:
        text (str): The text to search in
        pattern (str): The pattern to count
        d (int): Maximum mismatches allowed
        
    Returns:
        int: Count of approximate matches
    """
    count = 0
    k = len(pattern)
    
    for i in range(len(text) - k + 1):
        if hamming_distance(text[i:i + k], pattern) <= d:
            count += 1
    
    return count

def frequent_words_with_mismatches_and_reverse_complements(text, k, d):
    """
    Find most frequent k-mers considering both mismatches and reverse complements.
    
    For each k-mer Pattern, we count:
    Count_d(Text, Pattern) + Count_d(Text, Pattern_rc)
    
    Args:
        text (str): DNA string to analyze
        k (int): Length of k-mers
        d (int): Maximum mismatches allowed
        
    Returns:
        list: All k-mers with maximum combined count
    """
    freq_map = {}
    
    # Generate all possible k-mers
    all_kmers = set()
    for i in range(len(text) - k + 1):
        pattern = text[i:i + k]
        # Add both the pattern and all its neighbors to consider
        neighborhood = neighbors(pattern, d)
        all_kmers.update(neighborhood)
    
    # For each possible k-mer, calculate combined count
    for pattern in all_kmers:
        pattern_rc = reverse_complement(pattern)
        
        # Count approximate occurrences of pattern
        count_pattern = approximate_pattern_count(text, pattern, d)
        
        # Count approximate occurrences of reverse complement
        count_rc = approximate_pattern_count(text, pattern_rc, d)
        
        # Store the combined count
        freq_map[pattern] = count_pattern + count_rc
    
    # Find maximum combined count
    max_count = max(freq_map.values()) if freq_map else 0
    
    # Collect all k-mers with maximum count
    frequent_patterns = []
    for pattern in freq_map:
        if freq_map[pattern] == max_count:
            frequent_patterns.append(pattern)
    
    return frequent_patterns

# Read input from file
with open('data/part_2_data/sample_8.txt', 'r') as f:
    lines = f.read().strip().split('\n')
    text = lines[0]
    k, d = map(int, lines[1].split())

# Find frequent words with mismatches and reverse complements
result = frequent_words_with_mismatches_and_reverse_complements(text, k, d)

# Print results
print(' '.join(result))
