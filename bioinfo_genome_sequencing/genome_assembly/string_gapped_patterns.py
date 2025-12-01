def string_spelled_by_gapped_patterns(k, d, gapped_patterns):
    """
    Reconstruct a string from its sequence of (k,d)-mers.
    
    Args:
        k: Length of each k-mer in the pair
        d: Gap between the two k-mers
        gapped_patterns: List of (k,d)-mer pairs as tuples (first_kmer, second_kmer)
    
    Returns:
        Reconstructed string Text
    """
    # Split the gapped patterns into two separate lists
    first_patterns = [pattern[0] for pattern in gapped_patterns]
    second_patterns = [pattern[1] for pattern in gapped_patterns]
    
    # Reconstruct the prefix string from first patterns
    prefix_string = first_patterns[0]
    for i in range(1, len(first_patterns)):
        prefix_string += first_patterns[i][-1]
    
    # Reconstruct the suffix string from second patterns
    suffix_string = second_patterns[0]
    for i in range(1, len(second_patterns)):
        suffix_string += second_patterns[i][-1]
    
    # The full string is: prefix_string + gap + suffix_string
    # But they overlap, so we need to combine them correctly
    # prefix_string has length k + n - 1
    # suffix_string starts at position k + d
    # So we take the first (k + d) characters from prefix_string
    # and append suffix_string
    
    # Check if the overlapping parts match
    for i in range(len(second_patterns)):
        # Position in the final string where second pattern should start
        pos = k + d + i
        # Check if prefix_string[pos:pos+k] matches second_patterns[i]
        if pos + k <= len(prefix_string):
            if prefix_string[pos:pos+k] != second_patterns[i]:
                return None  # Patterns don't match
    
    # Construct the result
    result = prefix_string
    # Add any remaining characters from suffix_string
    overlap_start = k + d
    if overlap_start < len(prefix_string):
        # Verify the overlap matches
        overlap_length = len(prefix_string) - overlap_start
        if prefix_string[overlap_start:] != suffix_string[:overlap_length]:
            return None
        result = prefix_string + suffix_string[overlap_length:]
    else:
        # No overlap yet, need to add gap
        gap_size = overlap_start - len(prefix_string)
        result = prefix_string + suffix_string[0] * gap_size + suffix_string
    
    return result


def main():
    # Read input from file
    with open('bioinfo_genome_sequencing/datasets/dataset_10.txt', 'r') as f:
        lines = f.readlines()
    
    # Parse k and d from first line
    k, d = map(int, lines[0].strip().split())
    
    # Parse gapped patterns from second line
    pattern_strings = lines[1].strip().split()
    
    # Convert strings like "(GACC|GCGC)" to tuples ("GACC", "GCGC")
    gapped_patterns = []
    for pattern_str in pattern_strings:
        # Remove parentheses and split by |
        pattern_str = pattern_str.strip('()')
        parts = pattern_str.split('|')
        gapped_patterns.append((parts[0], parts[1]))
    
    # Reconstruct the string
    result = string_spelled_by_gapped_patterns(k, d, gapped_patterns)
    
    if result is None:
        print("Error: Cannot reconstruct string from given patterns")
        result = "Error: Invalid patterns"
    
    # Write output to file
    with open('bioinfo_genome_sequencing/datasets/output_10.txt', 'w') as f:
        f.write(result)
    
    print(f"Reconstructed string: {result}")
    print(f"Output written to output.txt")


if __name__ == "__main__":
    main()