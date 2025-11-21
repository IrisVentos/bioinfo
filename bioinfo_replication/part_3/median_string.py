#calculating the difference/distance between a pattern and one of the t DNA strings, and minimize that distance for all : Median String

from helpers import hamming_distance


def distance_between_pattern_and_strings(pattern, dna_strings):
    """
    Calculate the sum of distances between a pattern and each string in DNA collection.
    
    For each DNA string, finds the k-mer with minimum Hamming distance to the pattern,
    then sums these minimum distances across all DNA strings.
    
    Args:
        pattern: The pattern string to compare against
        dna_strings: List of DNA strings
    
    Returns:
        Total distance (sum of minimum distances from pattern to each DNA string)
    """
    k = len(pattern)
    total_distance = 0
    
    # For each DNA string
    for text in dna_strings:
        min_hamming_distance = float('inf')  # Initialize to infinity
        
        # Check every k-mer in this DNA string
        for i in range(len(text) - k + 1):
            kmer = text[i:i + k]
            
            # Calculate Hamming distance between pattern and this k-mer
            current_distance = hamming_distance(pattern, kmer)
            
            # Keep track of minimum distance for this DNA string
            if current_distance < min_hamming_distance:
                min_hamming_distance = current_distance
        
        # Add the minimum distance to total
        total_distance += min_hamming_distance
    
    return total_distance


def main():
    """
    Read input and calculate distance between pattern and DNA strings.
    """
    # Read input file
    with open('data/part_3_data/sample_6.txt', 'r') as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]
    
    # Parse input
    pattern = lines[0]
    dna_strings = lines[1].split()
    
    # Calculate distance
    result = distance_between_pattern_and_strings(pattern, dna_strings)
    
    # Output result
    print(result)


if __name__ == "__main__":
    main()