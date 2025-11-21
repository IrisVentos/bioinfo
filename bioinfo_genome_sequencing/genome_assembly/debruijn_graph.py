def debruijn_graph_from_string(k, text):
    """
    Construct the De Bruijn graph from a string.
    
    Args:
        k: The k-mer size
        text: The input genome string
    
    Returns:
        Dictionary representing adjacency list
    """
    graph = {}
    
    # Generate all k-mers from the text
    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]
        
        # Extract (k-1)-mer prefix and suffix
        prefix = kmer[:-1]  # First k-1 characters
        suffix = kmer[1:]   # Last k-1 characters
        
        # Add edge from prefix to suffix
        if prefix not in graph:
            graph[prefix] = []
        graph[prefix].append(suffix)
    
    return graph


def main():
    # Read input from file
    with open('datasets/dataset_4.txt', 'r') as f:
        lines = f.read().strip().split('\n')
    
    # Parse input
    k = int(lines[0])
    text = lines[1]
    
    # Build De Bruijn graph
    graph = debruijn_graph_from_string(k, text)
    
    # Build output
    output_lines = []
    for node in sorted(graph.keys()):
        neighbors = ' '.join(graph[node])
        line = f"{node}: {neighbors}"
        output_lines.append(line)
    
    # Write to output file
    with open('datasets/output_4.txt', 'w') as f:
        f.write('\n'.join(output_lines))
    
    # Also print to console
    for line in output_lines:
        print(line)


if __name__ == "__main__":
    main()