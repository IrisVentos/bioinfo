from helpers import debruijn_graph_from_kmers, eulerian_path, path_to_genome

def string_reconstruction(k, patterns):
    """
    Reconstruct a string from its k-mer composition.
    
    Args:
        k: Length of k-mers
        patterns: List of k-mers
    
    Returns:
        Reconstructed genome string
    """
    # Step 1: Construct De Bruijn graph from k-mers
    db_graph = debruijn_graph_from_kmers(patterns)
    
    # Step 2: Find Eulerian path in the graph
    path = eulerian_path(db_graph)
    
    # Step 3: Convert path to genome string
    text = path_to_genome(path)
    
    return text


def main():
    # Read input from file
    with open('bioinfo_genome_sequencing/datasets/dataset_8.txt', 'r') as f:
        lines = f.readlines()
    
    # Parse k value
    k = int(lines[0].strip())
    
    # Parse k-mers (they're on the second line, space-separated)
    patterns = lines[1].strip().split()
    
    # Reconstruct the genome
    result = string_reconstruction(k, patterns)
    
    # Write output to file
    with open('bioinfo_genome_sequencing/datasets/output_8.txt', 'w') as f:
        f.write(result)
    
    print(f"Reconstructed genome: {result}")
    print(f"Output written to output.txt")


if __name__ == "__main__":
    main()