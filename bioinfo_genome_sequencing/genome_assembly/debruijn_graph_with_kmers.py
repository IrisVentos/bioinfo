def debruijn_graph_from_kmers(kmers):
    """
    Construct the De Bruijn graph from a collection of k-mers.
    
    Args:
        kmers: List of k-mers (may contain duplicates)
    
    Returns:
        Dictionary representing adjacency list
    """
    graph = {}
    
    # For each k-mer, create an edge from prefix to suffix
    for kmer in kmers:
        prefix = kmer[:-1]  # First k-1 characters
        suffix = kmer[1:]   # Last k-1 characters
        
        # Add edge (duplicates are preserved)
        if prefix not in graph:
            graph[prefix] = []
        graph[prefix].append(suffix)
    
    return graph


def main():
    # Read input from file
    with open('datasets/dataset_5.txt', 'r') as f:
        content = f.read().strip()
    
    # Parse k-mers (space-separated)
    kmers = content.split()
    
    # Build De Bruijn graph
    graph = debruijn_graph_from_kmers(kmers)
    
    # Build output - sorted by node
    output_lines = []
    for node in sorted(graph.keys()):
        neighbors = ' '.join(graph[node])
        line = f"{node}: {neighbors}"
        output_lines.append(line)
    
    # Write to output file
    with open('datasets/output_5.txt', 'w') as f:
        f.write('\n'.join(output_lines))
    
    # Also print to console
    for line in output_lines:
        print(line)


if __name__ == "__main__":
    main()