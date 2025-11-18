def path_to_genome(kmers):
    """
    Reconstruct a string from its genome path.
    
    Args:
        kmers: List of k-mers where consecutive k-mers overlap by k-1 symbols
    
    Returns:
        Reconstructed genome string
    """
    if not kmers:
        return ""
    
    # Start with the first k-mer
    genome = kmers[0]
    
    # Add the last symbol of each subsequent k-mer
    for i in range(1, len(kmers)):
        genome += kmers[i][-1]
    
    return genome


def main():
    # Read input from file
    with open('datasets/dataset_2.txt', 'r') as f:
        content = f.read().strip()
    
    # Parse k-mers from input (space-separated)
    kmers = content.split()
    
    # Reconstruct genome
    result = path_to_genome(kmers)
    
    # Print result
    print(result)
    with open('datasets/output_2.txt', 'w') as f:
        f.write(result)

if __name__ == "__main__":
    main()
