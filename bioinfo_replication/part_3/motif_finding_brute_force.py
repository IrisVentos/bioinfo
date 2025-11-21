def hamming_distance(str1, str2):
    """
    Calculate the Hamming distance between two strings of equal length.
    Hamming distance = number of positions where characters differ.
    """
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))


def distance_pattern_to_string(pattern, dna_string):
    """
    Find the minimum Hamming distance between a pattern and any k-mer in a DNA string.
    This represents how "close" the pattern is to the best match in that string.
    """
    k = len(pattern)
    min_distance = float('inf')
    
    # Check every possible k-mer in the DNA string
    for i in range(len(dna_string) - k + 1):
        kmer = dna_string[i:i + k]
        distance = hamming_distance(pattern, kmer)
        min_distance = min(min_distance, distance)
    
    return min_distance


def distance_pattern_to_dna(pattern, dna_strings):
    """
    Calculate d(Pattern, Dna) = sum of distances from pattern to each DNA string.
    This is the total "cost" of choosing this pattern as our motif.
    """
    total_distance = 0
    for dna_string in dna_strings:
        total_distance += distance_pattern_to_string(pattern, dna_string)
    return total_distance


def generate_all_kmers(k):
    """
    Generate all possible k-mers using nucleotides A, C, G, T.
    For k=3: AAA, AAC, AAG, AAT, ACA, ACC, etc.
    Total: 4^k possible k-mers
    """
    nucleotides = ['A', 'C', 'G', 'T']
    
    if k == 1:
        return nucleotides
    
    # Recursively build k-mers by adding nucleotides to (k-1)-mers
    smaller_kmers = generate_all_kmers(k - 1)
    all_kmers = []
    
    for kmer in smaller_kmers:
        for nucleotide in nucleotides:
            all_kmers.append(kmer + nucleotide)
    
    return all_kmers


def median_string(k, dna_strings):
    """
    Find the k-mer that minimizes the total distance to all DNA strings.
    This is the "median string" - the most representative motif.
    """
    min_distance = float('inf')
    best_pattern = ""
    
    # Try every possible k-mer
    all_possible_kmers = generate_all_kmers(k)
    
    for pattern in all_possible_kmers:
        # Calculate how well this pattern matches all DNA strings
        current_distance = distance_pattern_to_dna(pattern, dna_strings)
        
        # Keep track of the best pattern so far
        if current_distance < min_distance:
            min_distance = current_distance
            best_pattern = pattern
    
    return best_pattern


def main():
    # Read input from file
    try:
        with open('data/part_3_data/sample_2.txt', 'r') as file:
            lines = file.read().strip().split('\n')
    except FileNotFoundError:
        print("Error: input.txt file not found!")
        print("Please create an input.txt file with the format:")
        print("Line 1: k (integer)")
        print("Line 2: space-separated DNA strings")
        return
    
    # Parse input
    k = int(lines[0])
    dna_strings = lines[1].split()
    
    print(f"Looking for motifs of length k = {k}")
    print(f"DNA strings: {dna_strings}")
    print(f"Total possible {k}-mers to check: 4^{k} = {4**k}")
    print()
    
    # Find the median string
    result = median_string(k, dna_strings)
    
    # Show detailed results
    print(f"Best motif found: {result}")
    
    # Show how well this motif matches each string
    print("\nDetailed analysis:")
    total_dist = 0
    for i, dna_string in enumerate(dna_strings):
        dist = distance_pattern_to_string(result, dna_string)
        total_dist += dist
        print(f"String {i+1} ({dna_string}): minimum distance = {dist}")
    
    print(f"\nTotal distance d({result}, Dna) = {total_dist}")


if __name__ == "__main__":
    main()